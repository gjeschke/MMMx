function [ecoor,atomtags,seq,options,correction,err] = RNA_loop_model(logfid,seq,HNP_lib,...
    anchor,acodes,transmat,ecodes,options,ntoffset,environ)
% function [ecoor,atomtags,seq,options,correction,err] = RNA_loop_model(logfid,seq,HNP_lib,...
%    anchor,acodes,transmat,ecodes,options,ntoffset,environ)
%
% Generates a model of an RNA loop with given sequence and anchor
% nucleotide(s) that does not clash with an environment
%
% Input:
%
% logfid        file identifier for log file
% seq           sequence in single letter code, allowed nucleotides A,C,G,U
%               modelled loop must be at least two nucleotides long
%               the code for the previous (anchor) nt and, if targeted,
%               next anchor nt must be included, use any base code for
%               non-existing anchors
% HNP_lib       HNP (Humphrey-Naranyan/Pyle) fragment library with fields
%               .fragments          full fragments
%               .shortfrag          short version of fragments
%               .non_clash_table    table of noc-clashing nucleotide pairs
% anchor        coordinates of the P(j),C4'(j),P(j+1) atoms of the
%               final anchor nucleotide, defaults to empty array, which is
%               interpreted as absence of a final anchor
% acodes        vector of fragment indices in the library that are allowed
%               for the initial anchor, defaults to all fragments in the
%               library; if not present or empty, no initial anchor is
%               assumed
% transmat      4x4 affine transformation matrix for rotation and
%               translation from the standard frame to the initial anchor
%               frame, defaults to identity transformation; if not present 
%               or empty, no initial anchor is assumed
% ecodes        vector of fragment indices in the library that are allowed
%               for the final nucleotide, used if the segment has a 
%               terminal anchor, defaults to all fragments in the library
% options       computation options, struct with fields
%               .attempts   maximum number of attempts for getting into the
%                           convergence radius of the final anchor,
%                           defaults to 100
%               .max_rot    maximum rotation of the loop end for matching 
%                           final anchor fragment, defaults to 8? per 4 nt
%               .anchor_acc RMSD accuracy of anchor coordinate fit,
%                           defaults to 0.25 ?
%               .fit_tol    contribution of a nucleotide to the convergence
%                           radius for a final anchor P atom, defaults to
%                           0.5 ?
%               .LoN        length per nt, defaults to 7 ?
%               .maxtime    maximum runtime in hours, defaults to 15 min
%               .anchor0    vector (1,4) defining an initial anchor of a 
%                           longer RNA segment by its nucleotide number and
%                           P atom coordinates; if present and if the loop
%                           is an insert, deformation of the whole segment
%                           is allowed for reaching the target anchor, the
%                           nucleotide number is relative to the first nt
%                           of the loop
%               .clash_thr  clash threshold, defaults to 1.4 Angstroem
% ntoffset      offset for the number of the first nucleotide, default: 0
% environ       coordinates of environment atoms with which the loop should
%               not clash, default to empty array
%
% Output:
%
% ecoor         extended coordinates (residue number, xyz) for atoms of the
%               RNA  model
% atomtags      cell vector of atom tag strings corresponding to ecoor
% seq           nucleotide sequence with the anchor nucleotides removed
% options       as the input options, but augmented by defaults
% correction    defines the transformation that superimposes the projected
%               C5'-terminal anchor fragment with the true C5'-terminal
%               anchor coordinates, structure with field
%               .trans      translation vector
%               .rot        vector defining a rotation axis (first 3 values)
%                           and rotation angle (in degree) about this axis
%               .length     length from first atom to anchor point
%               .Pinitial   coordinates of initial P atom  
% err           error code
%                 -1  clashes could not be corrected
%                 -2  run out of time
%                 -3  distance between anchors too large for number of
%                     nucleotides
%
% G. Jeschke, 2017-2021

atomtags = cell(0,1);
ecoor = [];
max_rpnt = 7; % maximum distance per nucleotide , 7 Angstroem
err = 0;
correction = [];

% Check whether the loop has two anchors, and if so, if the anchors are
% sufficiently close
if exist('transmat','var') && ~isempty(transmat) % there exists an initial anchor
    if exist('anchor','var') && ~isempty(anchor) % there exists a terminal anchor
        rpnt = norm(transmat(1:3,4)'-anchor(1,:))/(length(seq)-1); % distance per nucleotide
        if rpnt > max_rpnt
            seq = seq(2:end-1); % remove anchor nucleotides from sequence
            err = -3;
            return
        end
    end
end

if ~exist('environ','var')
    environ = [];
end

if ~exist('ntoffset','var') || isempty(ntoffset)
    ntoffset = 0;
end

if ~exist('options','var') || ~isfield(options,'clash_thr')
    options.clash_thr = 1.4;
end

if ~isfield(options,'maxtime') || isempty(options.maxtime)
    options.maxtime = 0.25;
end

if ~isfield(options,'attempts') || isempty(options.attempts)
    options.attempts = 100;
end

transmat0 = transmat;

tstart = tic;
runtime = toc(tstart);

while isempty(ecoor) && runtime < 3600*options.maxtime
    
    transmat = transmat0;
    
    [~,code,~,correction,lerr,statistics] = mk_RNA_loop_backbone(seq,...
        HNP_lib,anchor,acodes,transmat,ecodes,options);
    if lerr == 0
        % fprintf(logfid,'Successful RNA loop modelling after %5.1f s and %i trials\n',3600*statistics.runtime,statistics.trials);
    else
        fprintf(logfid,'RNA modelling failed after  %5.1f s and %i trials\n',3600*statistics.runtime,statistics.trials);
%         for ke = 1:length(statistics.errors)
%             if statistics.errors(ke) > 0.5/statistics.trials                
%                 fprintf(logfid,'For %6.2f%% of the attempts, %s\n',100*statistics.errors(ke),statistics.errmsg{ke});
%             end
%         end
        % the following messages apply only if there is a terminal anchor
        if statistics.min_convg < 1e5
            fprintf(logfid,'Best convergence was %4.2f %s\n',statistics.min_convg,char(197));
        end
        if statistics.amin_rmsd < 1e5
            fprintf(logfid,'Best anchor rmsd was %4.2f %s\n',statistics.amin_rmsd,char(197));
        end
        if statistics.amin_rot < 1e5
            fprintf(logfid,'Minimum rotation was %4.2f degree\n',statistics.amin_rot);
        end
        if statistics.amin_shift < 1e5
            fprintf(logfid,'Minimum translation was %4.2f %s\n',statistics.amin_shift,char(197));
        end
        runtime = toc(tstart);
        % if backbone generation was unsuccessful, but there is still time,
        % try again
        continue
    end
    
    if ~exist('transmat','var') || isempty(transmat)
        transmat = eye(4);
    end

    [ecoor,atomtags] = mk_RNA(HNP_lib.fragments,seq(1:end-1),code(1:end-1),HNP_lib.non_clash_table,ntoffset,transmat);
       
    if isempty(ecoor)
        % fprintf(logfid,'RNA backbone could not be decorated with bases.\n');
        runtime = toc(tstart);
        % if full RNA generation was unsuccessful, but there is still time,
        % try again
        continue
    else
        valid = find(ecoor(:,1)>ecoor(1,1)); % atoms not in the initial anchor nt
        atomtags = atomtags(valid);
        ecoor = ecoor(valid,:);
    end

    % distribute correction of translation and rotation of the final
    % nucleotide (if any) over the whole loop
    if correction.corrected && ~isfield(options,'anchor0')
        ecoor = anchor_correction(ecoor,correction);
    end

    % try to repair clashes, if requested
    if ~isempty(environ) % && ~isfield(options,'anchor0')
        coor = ecoor(:,2:4);
        if ~isempty(anchor)
            fixed = anchor(1,:);
        else
            fixed = zeros(0,3);
        end
        if ~isempty(transmat)
            f2 = transmat(1:3,4);
            fixed = [f2'; fixed]; %#ok<AGROW>
        end
        options.max_shift = 2;
        [coor,min_dist,~,max_shift] = clash_repair(coor,environ,fixed,options);
        ecoor(:,2:4) = coor;
        if min_dist < options.clash_thr || max_shift > options.max_shift
            % fprintf(logfid,'Clash repair failed\n');
            ecoor = [];
            err = -1;
        else
            fprintf(logfid,'Clash repair led to a maximum shift of %6.2f %s\n',max_shift,char(197));
            err = 0;
        end
    end
    runtime = toc(tstart);
    
end

if isempty(ecoor) && runtime >= 3600*options.maxtime
    err = -2;
end

seq = seq(2:end-1);



function [coor,code,transmat,correction,err,statistics] = mk_RNA_loop_backbone(seq,HNP_lib,anchor,acodes,transmat,ecodes,options)
% function [coor,code,transmat,correction,err,statistics] = mk_RNA_loop_backbone(seq,HNP_lib,anchor,acodes,transmat,ecodes,options)
%
% Generates an RNA backbone by random selection of conformers from a 
% fragment library, 
% see also E. Humphris-Narayanan, A. M. Pyle, J. Mol. Biol. 2012, 421, 6-26
% the segment attaches to an initial anchor and leads to a final anchor
% either of the anchors or both may be missing
%
% seq           sequence in single letter code, allowed nucleotides A,C,G,U
%               modelled loop must be at least two nucleotides long
%               the code for the previous (anchor) nt and, if targeted,
%               next anchor nt must be included, use any base code for
%               non-existing anchors
% HNP_lib       HNP (Humphrey-Naranyan/Pyle) fragment library with fields
%               .fragments          full fragments
%               .shortfrag          short version of fragments
%               .non_clash_table    table of noc-clashing nucleotide pairs
% anchor        coordinates of the P(j),C4'(j),P(j+1) atoms of the
%               final anchor nucleotide, defaults to empty array, which is
%               interpreted as absence of a final anchor
% acodes        vector of fragment indices in the library that are allowed
%               for the initial anchor, defaults to all fragments in the
%               library, if not present or empty, no initial anchor is
%               assumed
% transmat      4x4 affine transformation matrix for rotation and
%               translation from the standard frame to the initial anchor
%               frame, defaults to identity transformation, if not present 
%               or empty, no initial anchor is assumed
% ecodes        vector of fragment indices in the library that are allowed
%               for the final nucleotide, used if the segment attaches to a
%               stem, optional, defaults to all fragments in the library
% options       computation options, struct with fields
%               .attempts   maximum number of attempts for getting into the
%                           convergence radius of the final anchor,
%                           defaults to 100
%               .max_rot    maximum rotation of the loop end for matching 
%                           final anchor fragment, defaults to 20? per 4 nt
%               .anchor_acc RMSD accuracy of anchor coordinate fit,
%                           defaults to 0.25 ?
%               .fit_tol    contribution of a nucleotide to the convergence
%                           radius for a final anchor P atom, defaults to
%                           0.2 ?
%               .LoN        length per nt, defaults to 7 ?
%               .maxtime    maximum runtime in hours, defaults to 25 min
%               .anchor0    vector (1,4) defining an initial anchor of a 
%                           longer RNA segment by its nucleotide number and
%                           P atom coordinates; if present and if the loop
%                           is an insert, deformation of the whole segment
%                           is allowed for reaching the target anchor, the
%                           nucleotide number is relative to the first nt
%                           of the loop
%
% coor          [2*N+2,3] matrix of coordinates of atoms, where N is the
%               length of seq, the additional atom at the end is the P atom
%               of a hypothetical next nt, the additional atom at the
%               beginning the C4' atom of the previous nt
% code          vector of fragment numbers in the library, length(seq),
%               since the initial and final anchor residues are assigned a 
%               code, if no target anchor exists, the last code is zero
% transmat      matrix for affine transformation of a model in the standard
%               frame to the frame defined by the C3'-terminal anchor 
%               nucleotide
% correction    defines the transformation that superimposes the projected
%               C5'-terminal anchor fragment with the true C5'-terminal
%               anchor coordinates, structure with field
%               .trans      translation vector
%               .rot        vector defining a rotation axis (first 3 values)
%                           and rotation angle (in degree) about this axis
%               .length     length from first atom to anchor point
%               .Pinitial   coordinates of initial P atom
% err           error code if no model could be produced
%               positive numbers: residue number, where the chain did not
%               find back to the convergence radius
%               -1  P atom of anchor residue was not reached
%               -2  orientation of anchor residue could not be matched
%               -3  no final nt that is both non-clashing and allowed for
%                   target 
% statistics    struct with fields
%               .runtime    runtime in hours
%               .trials     number of trials used
%               .errnum     axis for error code
%               .errors     error code statistics
%               .errmsg     error messages
%               
%
% G. Jeschke, 1/2017-25.12.2017

% initialize error messages

statistics.errmsg{1} = 'P atom of terminal anchor was not reached';
statistics.errmsg{2} = 'orientation of terminal anchor could not be matched';
statistics.errmsg{3} = 'no final nucleotide was allowed and non-clashing';

verbose = false;

correction.trans = [0,0,0]; % by default, no correction translation
correction.rot = eye(3); % by default, no correction rotation
correction.length = 0;
correction.corrected = false;

if ~isempty(transmat)
    initial_anchored = true;
else
    initial_anchored = false;
end

slc = 'ACGU'; 

ntinit = seq(1); % nucleotide code for C3' (previous segment) anchor
nttarget = seq(end); % nucleotide code for C5' (next segment) anchor
seq = seq(1:end-1); % sequence for only the modelled loop segment

if ~isfield(HNP_lib,'non_clash_table') || isempty(HNP_lib.non_clash_table)
    clash_test = false;
else
    clash_test = true;
end

if ~exist('anchor','var')
    anchor = [];
end

if ~exist('acodes','var') || isempty(acodes)
    acodes = 1:length(HNP_lib.shortfrag);
    initial_anchored = false;
end

if ~exist('transmat','var') || isempty(transmat)
    transmat = eye(4);
    initial_anchored = false;
end

if ~exist('ecodes','var') || isempty(ecodes)
    ecodes = 1:length(HNP_lib.shortfrag);
end

if ~exist('options','var') || ~isfield(options,'attempts')
    options.attempts = 100; % number of attempts before declaring failure
end

if ~exist('options','var') || ~isfield(options,'max_rot')
    options.max_rot = 10; % maximum rotation of the loop end for matching final anchor fragment (per 4 nt)
end

if ~exist('options','var') || ~isfield(options,'anchor_acc')
    options.anchor_acc = 0.5; % RMSD accuracy of anchor coordinate fit
end

if ~exist('options','var') || ~isfield(options,'fit_tol')
    options.fit_tol = 0.5; % tolerance per nucleotide for convergence radius of target P atom
end

if ~isfield(options,'maxtime')
    options.maxtime = 0.25; % maximum runtime in hours
end

if ~isfield(options,'LoN')
    options.LoN = 7; % maximum length per nucleotide in ?
end

if isfield(options,'anchor0')
    len = length(seq) - options.anchor0(1) - 2;
else
    len = length(seq);
end
fit_thresh = options.fit_tol*(len+2);
if fit_thresh < 1
    fit_thresh = 1;
end
max_rot = len*options.max_rot/4;
if max_rot < 10
    max_rot = 10;
end

error_statistics = zeros(1,4+length(seq));
statistics.errnum = -3:length(seq);
err = -100;
trials = 0;
tstart = tic;
runtime = toc(tstart);
if length(seq) < 4
    maxtrials = 4*length(acodes)*length(HNP_lib.shortfrag)^(length(seq)-2);
else
    maxtrials = 1e6;
end

if ~isempty(anchor)
    Ptarget = anchor(1,:); % target coordinate of P atom
    if initial_anchored
        Pinitial = transmat(1:3,4)';
        % determine if loop ends at a final anchor and, if so, store target
        % coordinate
        rP5p = norm(Pinitial-Ptarget); % distance to the target P atom
        det = rP5p - length(seq)*options.LoN; % P atom must within contour length to the target
        if det > fit_thresh
            err = -1;
            statistics.errors(3) = 1;
            return
        end
    end
end

bcode = zeros(1,length(seq)+1);
bcode(1) = strfind(slc,upper(ntinit));
for k = 2:length(seq)-1
    base = upper(seq(k));
    bcode(k) = strfind(slc,base);
end

statistics.min_convg = 1e6;
statistics.amin_rmsd = 1e6;
statistics.amin_rot = 1e6;
statistics.amin_shift = 1e6;
min_approaches = 1e6*ones(1,length(seq)-1);

in_clash_test_range = initial_anchored;

while err ~= 0 && runtime < 3600*options.maxtime && trials < maxtrials

    nl = length(HNP_lib.shortfrag);
    coor = zeros(2*length(seq)+2,3);
    code = zeros(1,length(seq)+1);
    
    % select initial fragment from the ones that are allowed and get its
    % standard coordinates
    na = length(acodes);
    acode = 1 + floor(na*rand-eps);
    code(1) = acodes(acode);
    sfrag = [HNP_lib.fragments(code(1)).(ntinit).coor(HNP_lib.fragments(code(1)).(ntinit).assign.previous(2:3),:);...
        HNP_lib.fragments(code(1)).(ntinit).coor(HNP_lib.fragments(code(1)).(ntinit).assign.next(2:3),:);];
    
    % amend to coordinate array for affine transformation and transform to
    % initial anchor frame
    sfrag = [sfrag ones(4,1)]; %#ok<AGROW>
    sfrag = transmat*sfrag';
    sfrag = sfrag';
    
    % store pseudo-torsion angle defining coordinates of initial nucleotide
    coor(1:4,:) = sfrag(:,1:3);
    
    
    % store coordinate of initial P atom, in the anchor0 case, the O3' atom
    % is used
    Pinitial = coor(1,:);
    if isfield(options,'anchor0')
        correction.Pinitial = options.anchor0(2:4);
    else
        correction.Pinitial = Pinitial;
    end
        
    % determine if loop ends at a final anchor and, if so, store target
    % coordinate
    if ~isempty(anchor)
        Ptarget = anchor(1,:); % target coordinate of P atom
        targeted = true;
        rP5p = norm(Pinitial-Ptarget); % distance to the target P atom
        det = rP5p - length(seq)*options.LoN; % P atom must within contour length to the target
        % if the target cannot be reached, return empty coordinates
        if initial_anchored && det > fit_thresh  
            err = -1;
            code = [];
            coor = [];
            runtime = toc(tstart);
            trials = trials+1;
            error_statistics(4+err) = error_statistics(4+err) + 1;
            continue;
        end
    else
        targeted = false;
    end
    
    poi = 2; % the first two atom coordinates will not change
    % for all except the last inserted nucleotide, use a Monte Carlo procedure
    passed = true;
    
    for k = 2:length(seq)-1
        det = 1e6; % initialize distance from target
        count = 0;
        coor1 = coor(poi:poi+2,:); % atom coordinates defining local frame
        [Rp,orig] = get_trafo(coor1); % local frame
        origmat = repmat(orig,4,1);
        % try to make a step that brings chain into current convergence radius
        while det > fit_thresh && count < options.attempts
            if in_clash_test_range && clash_test % allow only for non-clashing pair with previous nt
                pair = 4*(bcode(k-1)-1) + bcode(k);
                non_clash_table = HNP_lib.non_clash_table{pair};
                non_clashing = non_clash_table{code(k-1)};
                code_pointer = 1 + floor(length(non_clashing)*rand-eps);
                code(k) = non_clashing(code_pointer);
            else % allow for any new fragment
                code(k) = 1 + floor(nl*rand-eps);
            end
            coor2 = HNP_lib.shortfrag{code(k)}; % coordinates of current fragment
            coor2 = coor2*Rp + origmat; % transform to local frame
            if targeted % check whether inside convergence radius
                if initial_anchored
                    rP5p = norm(coor2(3,:)-Ptarget); % distance to the target P atom
                    det = rP5p - (length(seq)-k)*options.LoN; % P atom must within contour length to the target
                else
                    det = 0;
                end
            else
                det = 0;
            end
            count = count + 1;
        end
        if det - fit_thresh < min_approaches(k)
            min_approaches(k) = det - fit_thresh;
        end
        % if the target was not reached, return empty coordinates
        if det > fit_thresh
            err = k;
            code = [];
            coor = [];
            runtime = toc(tstart);
            trials = trials+1;
            passed = false;
            error_statistics(4+err) = error_statistics(4+err) + 1;
            break;
        end
        coor(poi+1:poi+4,:) = coor2; % store current fragment coordinates
        poi = poi + 2; % count up coordinate pointer
        in_clash_test_range = true; % switch on clash tests if requested
    end
    % if something went wrong and there is still time, make another attempt
    if ~passed
        continue;
    end
    
    coor1 = coor(poi:poi+2,:); % cut off coordinate array
    
    k = length(seq);
    
    if targeted
        correction.length = norm(Ptarget-Pinitial);
        correction.corrected = true;
        [Rp,orig] = get_trafo(coor1);
        % for the last nucleotide, look for the best fit of all fragments
        base = upper(seq(k));
        available = ecodes;
        if clash_test
            bcode(k) = strfind(slc,base);
            pair = 4*(bcode(k-1)-1) + bcode(k);
            non_clash_table = HNP_lib.non_clash_table{pair};
            available = non_clash_table{code(k-1)};
            apoi = 0;
            % check which non-clashing fragments are allowed to target
            for ka = 1:length(available)
                if min(abs(ecodes-available(ka))) < eps
                    apoi = apoi + 1;
                    available(apoi) = available(ka);
                end
            end
            available = available(1:apoi);
            if isempty(available) % no final fragment that is both non-clashing and allowed to target
                err = -3;
                code = [];
                coor = [];
                runtime = toc(tstart);
                trials = trials+1;
                error_statistics(4+err) = error_statistics(4+err) + 1;
                continue;
            end
        end
        rP5_vec = 1e6*ones(1,length(HNP_lib.fragments));
        for kf = available
            coor2 = HNP_lib.shortfrag{kf};
            coor2 = coor2*Rp + repmat(orig,4,1);
            rP5_vec(kf) = norm(coor2(3,:)-Ptarget); % distance to the target P atom
        end
        [mi,bf] = min(rP5_vec);
        
        err = length(seq);
        if mi < statistics.min_convg
            statistics.min_convg = mi;
        end
        if mi < fit_thresh || ~initial_anchored % P atom within convergence radius
            coor2 = HNP_lib.shortfrag{bf};
            coor2 = coor2*Rp + repmat(orig,4,1);
            code(k) = bf;
            coor(poi+1:poi+4,:) = coor2;
            % find best matching anchor fragment (minimum rotation at target)
            available = 1:length(HNP_lib.fragments);
            if clash_test
                % which fragments are non-clashing
                bc5p = strfind(slc,nttarget);
                pair = 4*(bcode(k)-1) +bc5p;
                non_clash_table = HNP_lib.non_clash_table{pair};
                available = non_clash_table{code(k-1)};
            end
            coort = coor2(2:4,:);
            min_rot = max_rot + eps;
            for kf = available
                fanchor = HNP_lib.fragments(kf).A.coor([HNP_lib.fragments(kf).A.assign.previous(2:3) HNP_lib.fragments(kf).A.assign.next(2)],:);
                [~, ~, transmat0] = superimpose_3points(anchor,fanchor);
                linkage = HNP_lib.fragments(kf).A.coor(HNP_lib.fragments(kf).A.assign.previous,:);
                linkage = [linkage ones(3,1)]*transmat0';
                [rmsd, ~, transmat_link] = superimpose_3points(linkage(:,1:3),coort);
                EV = affine2EV(transmat_link);
                fail = rmsd > 2*options.anchor_acc;
                fail = fail + (norm(transmat_link(1:3,4)) > fit_thresh);
                fail = fail + (abs(EV(4)) > max_rot);
                if verbose && fail <= 1
                    if rmsd > 2*options.anchor_acc
                        fprintf(2,'RMSD: %4.2f ?, ',rmsd);
                    else
                        fprintf(1,'RMSD: %4.2f ?, ',rmsd);
                    end
                    if norm(transmat_link(1:3,4)) > fit_thresh
                        fprintf(2,'translation %4.2f ?, ',norm(transmat_link(1:3,4)));
                    else
                        fprintf(1,'translation %4.2f ?, ',norm(transmat_link(1:3,4)));
                    end
                    if abs(EV(4)) > max_rot
                        fprintf(2,'rotation %4.1f?\n',EV(4));
                    else
                        fprintf(1,'rotation %4.1f?\n',EV(4));
                    end
                end
                if ~initial_anchored || ...
                        (abs(EV(4)) < min_rot && rmsd < 2*options.anchor_acc && norm(transmat_link(1:3,4)) < fit_thresh)
                    err = 0;
                    correction.trans = transmat_link(1:3,4)';
                    correction.rot = EV;
                    min_rot = abs(EV(4));
                    code(k+1) = kf;
                    if initial_anchored
                        correction.local = true;
                    else
                        correction.local = false;
                    end
                else
                    if abs(EV(4)) < statistics.amin_rot
                        statistics.amin_rot = abs(EV(4));
                    end
                    if rmsd < statistics.amin_rmsd
                        statistics.amin_rmsd = rmsd;
                    end
                    if norm(transmat_link(1:3,4)) < statistics.amin_shift
                        statistics.amin_shift = norm(transmat_link(1:3,4));
                    end
                end
            end
        else
            err = -1;
            code = [];
            coor = [];
        end
    else % non-targeted
        [Rp,orig] = get_trafo(coor1);
        % for the last nucleotide, look for the best fit of all fragments
        base = upper(seq(k));
        available = ecodes;
        if clash_test
            bcode(k) = strfind(slc,base);
            pair = 4*(bcode(k-1)-1) + bcode(k);
            non_clash_table = HNP_lib.non_clash_table{pair};
            available = non_clash_table{code(k-1)};
            apoi = 0;
            % check which non-clashing fragments are allowed to target
            for ka = 1:length(available)
                if min(abs(ecodes-available(ka))) < eps
                    apoi = apoi + 1;
                    available(apoi) = available(ka);
                end
            end
            available = available(1:apoi);
            if isempty(available) % no final fragment that is both non-clashing and allowed to target
                err = -3;
                code = [];
                coor = [];
                runtime = toc(tstart);
                trials = trials+1;
                error_statistics(4+err) = error_statistics(4+err) + 1;
                continue;
            end
        end
        code_pointer = 1 + floor(length(available)*rand-eps);
        code(k) = available(code_pointer);
        coor2 = HNP_lib.shortfrag{code(k)}; % coordinates of current fragment
        coor2 = coor2*Rp + repmat(orig,4,1); % transform to local frame
        coor(poi+1:poi+4,:) = coor2; % store current fragment coordinates
        Pfinal = coor2(1,:);
        correction.length = norm(Pfinal-Pinitial);
        err = 0;
    end
    runtime = toc(tstart);
    trials = trials+1;
    error_statistics(4+err) = error_statistics(4+err) + 1;
end

statistics.runtime = runtime/3600;
statistics.trials = trials;
statistics.errors = error_statistics/trials;

function [Rp,orig] = get_trafo(coor)

orig = coor(2,:);
coor = coor - repmat(orig,3,1);
x = coor(1,:)-coor(2,:); 
x = x/norm(x);    % unit vector along x
yp = coor(3,:)-coor(2,:); 
yp = yp/norm(yp);
z = cross_rowvec(x,yp); % z axis is perpendicular on xy plane
z = z/norm(z);
y = cross_rowvec(z,x); % real (corrected) y axis
Rp = [x;y;z];

function c=cross_rowvec(a,b)
% A fast cross product that works only for two three-element row vectors 

c = [a(2)*b(3)-a(3)*b(2),a(3)*b(1)-a(1)*b(3),a(1)*b(2)-a(2)*b(1)];

function ecoor = anchor_correction(ecoor,corr)
% function ecoor = anchor_correction(ecoor,corr)
% 
% Correction of an extended coordinate set of an RNA loop for matching the
% final anchor, use mk_RNA_loop_backbone.m and mk_RNA.m to obtain
% corr and ecoor
% the correcting translation and rotation is proportional to the distance
% between the current atom and the initial P atom, divided by the distance
% between the target P atom and the initial P atom
% this distributes the deviation between target and final P atom of the 
% original construct, so that each individual bond length, bond angle, and 
% torsion angle is only weakly affected
%
% ecoor [m,4] coordinate array for m atoms, first column is residue
%       number, columns 2-4 are Cartesian coordinates
% corr  correction parameters, structure   
%       .trans      translation vector
%       .rot        4-element vector, elements 1:3 specify rotation axis,
%                   element 4 specifies rotation angle in degrees
%       .length     length from first atom to anchor point
%       .Pinitial   coordinates of initial P atom      
%       .local      Boolean, local correction if true, otherwise global
%                   correction
%
% G. Jeschke, 26.12.2017

[m,~] = size(ecoor);
for k = 1:m
    if corr.local
        progress = norm(ecoor(k,2:4)-corr.Pinitial)/corr.length;
    else
        progress = 1;
    end
    tm1 = affine('translation',progress*corr.trans);
    tm2 = affine('rotn',progress*pi*corr.rot(4)/180,corr.rot(1:3));
    tm = tm1*tm2;
    ccoor = [ecoor(k,2:4) 1]*tm';
    ecoor(k,2:4) = ccoor(1:3);
end
