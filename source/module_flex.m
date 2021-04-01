function [entity,exceptions,failed] = module_flex(control,logfid,entity)
%
% MMModel Run modelling pipeline
%
%   [exceptions,entity] = MMModel(control_file)
%   Parses control file for modelling directives and arguments and runs the
%   specified modelling pipeline
%
% INPUT
% control       control structure with fields
%               .name           'ensemblefit', unused
%               .options        options struct
%               .directives     array of struct with directives
%               .entity_output  number of the entity to be output
% logfid        file handle for log file
% entity        input entity, can be empty for free-standing flexible
%               domains, defaults to empty
%
% OUTPUT
% entity        output entity, empty in NOFIT mode
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
% failed        flag indicated whther the module failed
%
% SUPPORTED DIRECTIVES
%
% parallel      number of parallel runs per block
% interactive   if present, information on progress is displayed
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
exceptions = {[]};
warnings = 0;
failed = false;

if ~exist('entity','var')
    entity = [];
end
% set defaults

keep_sidegroup_clashes = false; % delete models where side groups clash

opt.parnum = 100; % number of trials performed in the parfor loop
opt.disp_update = 200; % cycles between display updates in interactive mode
opt.interactive = false;
maxlen = 1000; % maximum expected loop length for memory pre-allocation
M_max = 10; % maximum factor for rejection sampling (default)
f_update = 0.001; % update factor for sampling functions in rejection sampling (default)

fname = 'mmmx_flex';
pdbid = 'MMMX';
chainid = 'A';

restraints.ddr(1).labels{1} = '';
restraints.oligomer(1).label = '';
restraints.depth(1).label = '';
restraints.aprop = zeros(1,maxlen);
restraints.bprop = zeros(1,maxlen);
restraints.cprop = zeros(1,maxlen);

ddr_poi = 0;
oligomer_poi = 0;
depth_poi = 0;
rejection_sampling = false;
% read restraints
for d = 1:length(control.directives)
    switch lower(control.directives(d).name)
        case 'sequence'
            % first, last residue number and sequence
            restraints.initial =  str2double(control.directives(d).options{1});
            restraints.final =  str2double(control.directives(d).options{2});
            restraints.sequence = control.directives(d).options{3};
        case 'interactive'
            opt.interactive = true;
            if ~isempty(control.directives(d).options) && ~isempty(control.directives(d).options{1})
                opt.disp_update = str2double(control.directives(d).options{1});
            end
        case 'loose'
            keep_sidegroup_clashes = true;
        case 'parallel'
            opt.parnum = str2double(control.directives(d).options{1});
        case 'save'
            fname = control.directives(d).options{1};
            if length(control.directives(d).options) >= 2
                pdbid = control.directives(d).options{2};
                if length(pdbid) > 4
                    pdbid = pdbid(1:4);
                elseif length(pdbid) < 4
                    pdbid = pad(pdbid,'*');
                end
            end
            if length(control.directives(d).options) >= 3
                chainid = control.directives(d).options{3};
                if length(chainid) > 1
                    chainid = chainid(1);
                end
            end
        case 'c_anchor'
            restraints.c_anchor = control.directives(d).options{1};
        case 'n_anchor'
            restraints.n_anchor = control.directives(d).options{1};
        case 'a_prop'
            [nr,~] = size(control.directives(d).block);
            if nr > 0 % at least one restraint in this block
                for kr = 1:nr
                    res = str2double(control.directives(d).block{kr,1});
                    prop = str2double(control.directives(d).block{kr,2});
                    restraints.aprop(res) = prop;
                end
            else % empty block
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_flex:empty_a_prop_block', 'a_prop block %i is empty',d);
                record_exception(exceptions{warnings},logfid);
            end
        case 'b_prop'
            [nr,~] = size(control.directives(d).block);
            if nr > 0 % at least one restraint in this block
                for kr = 1:nr
                    res = str2double(control.directives(d).block{kr,1});
                    prop = str2double(control.directives(d).block{kr,2});
                    restraints.aprop(res) = prop;
                end
            else % empty block
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_flex:empty_b_prop_block', 'b_prop block %i is empty',d);
                record_exception(exceptions{warnings},logfid);
            end
        case 'c_prop'
            [nr,~] = size(control.directives(d).block);
            if nr > 0 % at least one restraint in this block
                for kr = 1:nr
                    res = str2double(control.directives(d).block{kr,1});
                    prop = str2double(control.directives(d).block{kr,2});
                    restraints.cprop(res) = prop;
                end
            else % empty block
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_flex:empty_c_prop_block', 'c_prop block %i is empty',d);
                record_exception(exceptions{warnings},logfid);
            end
        case 'ddr'
            ddr_poi = ddr_poi + 1; % increase ddr block counter
            fprintf(logfid,'ddr %s',control.directives(d).options{1}); % echo directive to log file
            restraints.ddr(ddr_poi).labels{1} = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % different labels
                restraints.ddr(ddr_poi).labels{2} = control.directives(d).options{2};
                fprintf(logfid,' %s\n',control.directives(d).options{2});
            else % same label at both sites
                restraints.ddr(ddr_poi).labels{2} = control.directives(d).options{1};
                fprintf(logfid,'\n');
            end
            if length(control.directives(d).options) > 2 % specification of maximum M for rejection sampling
                restraints.ddr(ddr_poi).M_max = str2double(control.directives(d).options{3});
            else
                restraints.ddr(ddr_poi).M_max = M_max;
            end
            if length(control.directives(d).options) > 3 % specification of sampling update factor
                restraints.ddr(ddr_poi).f_update = str2double(control.directives(d).options{4});
            else
                restraints.ddr(ddr_poi).f_update = f_update;
            end
            [nr,args] = size(control.directives(d).block);
            if nr > 0 % at least one restraint in this block
                restraints.ddr(ddr_poi).site1{nr} = ' ';
                restraints.ddr(ddr_poi).site2{nr} = ' ';
                restraints.ddr(ddr_poi).r(nr) = 0;
                restraints.ddr(ddr_poi).sigr(nr) = 0;
                restraints.ddr(ddr_poi).file{nr} = '*';
            else % empty block
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_flex:empty_ddr_block', 'ddr block %i is empty',d);
                record_exception(exceptions{warnings},logfid);
            end
            if args < 3 % misformed restraint line
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_flex:misformed_ddr', 'ddr restraint has less than three arguments');
                record_exception(exceptions{warnings},logfid);
                failed = true;
                return
            end
            for kr = 1:nr
                restraints.ddr(ddr_poi).site1{kr} = control.directives(d).block{kr,1};
                restraints.ddr(ddr_poi).site2{kr} = control.directives(d).block{kr,2};
                restraints.ddr(ddr_poi).file{kr} = '';
                arg3 = control.directives(d).block{kr,3};
                if arg3(1) == '@'
                    restraints.ddr(ddr_poi).file(kr) = arg3(2:end);
                    restraints.ddr(ddr_poi).r(kr) = [];
                    restraints.ddr(ddr_poi).sigr(kr) = [];
                    rejection_sampling = true;
                else
                    restraints.ddr(ddr_poi).r(kr) = str2double(arg3);
                    restraints.ddr(ddr_poi).sigr(kr) = str2double(control.directives(d).block{kr,4});
                end
                if args > 4
                    arg5 = control.directives(d).block{kr,5};
                    if ~isempty(arg5) && arg5(1) == '@'
                        restraints.ddr(ddr_poi).file{kr} = arg5(2:end);
                        rejection_sampling = true;
                    end
                end
                for karg = 1:args
                    fprintf(logfid,'  %s',control.directives(d).block{kr,karg});
                end
                fprintf(logfid,'\n');
            end
            fprintf(logfid,'\n\n');
        case 'oligomer'
            oligomer_poi = oligomer_poi + 1; % increase oligomer block counter
            fprintf(logfid,'oligomer %s',control.directives(d).options{1}); % echo directive to log file
            restraints.oligomer(oligomer_poi).label = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % multiplicity provided
                restraints.oligomer(oligomer_poi).multiplicity = str2double(control.directives(d).options{2});
                fprintf(logfid,' %s\n',control.directives(d).options{2});
            else % homodimer assumed
                restraints.oligomer(oligomer_poi).multiplicity = 2;
                fprintf(logfid,'\n');
            end
            if length(control.directives(d).options) > 2 % maximum M for rejection sampling provided
                restraints.oligomer(oligomer_poi).M_max = str2double(control.directives(d).options{3});
            else
                restraints.oligomer(oligomer_poi).M_max = M_max;
            end
            if length(control.directives(d).options) > 3 % maximum M for rejection sampling provided
                restraints.oligomer(oligomer_poi).f_update = str2double(control.directives(d).options{4});
            else
                restraints.oligomer(oligomer_poi).f_update = f_update;
            end
            [nr,args] = size(control.directives(d).block);
            if nr > 0 % at least one restraint in this block
                restraints.oligomer(oligomer_poi).site{nr} = ' ';
                restraints.oligomer(oligomer_poi).r(nr) = 0;
                restraints.oligomer(oligomer_poi).sigr(nr) = 0;
                restraints.oligomer(oligomer_poi).file{nr} = '*';
            else % empty block
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_flex:empty_oligomer_block', 'oligomer block %i is empty',d);
                record_exception(exceptions{warnings},logfid);
            end
            if args < 3 % misformed restraint line
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_flex:misformed_oligomer_restraint', 'oligomer restraint has less than three arguments');
                record_exception(exceptions{warnings},logfid);
                failed = true;
                return
            end
            for kr = 1:nr
                restraints.oligomer(oligomer_poi).site{kr} = control.directives(d).block{kr,1};
                restraints.oligomer(oligomer_poi).file{kr} = '';
                arg2 = control.directives(d).block{kr,2};
                if arg2(1) == '@'
                    restraints.oligomer(oligomer_poi).file(kr) = arg2(2:end);
                    restraints.oligomer(oligomer_poi).r(kr) = [];
                    restraints.oligomer(oligomer_poi).sigr(kr) = [];
                    rejection_sampling = true;
                else
                    restraints.oligomer(oligomer_poi).r(kr) = str2double(arg2);
                    restraints.oligomer(oligomer_poi).sigr(kr) = str2double(control.directives(d).block{kr,3});
                end
                if args > 3
                    arg4 = control.directives(d).block{kr,4};
                    if ~iempty(arg4) && arg4(1) == '@'
                        restraints.oligomer(oligomer_poi).file{kr} = arg4(2:end);
                        rejection_sampling = true;
                    end
                end
                for karg = 1:args
                    fprintf(logfid,'  %s',control.directives(d).block{kr,karg});
                end
                fprintf(logfid,'\n');
            end
            fprintf(logfid,'\n\n');
        case 'depth'
            depth_poi = depth_poi + 1; % increase depth block counter
            fprintf(logfid,'depth %s',control.directives(d).options{1}); % echo directive to log file
            restraints.depth(depth_poi).label = control.directives(d).options{1};
            [nr,args] = size(control.directives(d).block);
            if nr > 0 % at least one restraint in this block
                restraints.depth(depth_poi).site{nr} = ' ';
                restraints.depth(depth_poi).r(nr) = 0;
                restraints.depth(depth_poi).sigr(nr) = 0;
            else % empty block
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_flex:empty_depth_block', 'depth block %i is empty',d);
                record_exception(exceptions{warnings},logfid);
            end
            if args < 3 % misformed restraint line
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_flex:misformed_depth_restraint', 'depth restraint has less than three arguments');
                record_exception(exceptions{warnings},logfid);
                failed = true;
                return
            end
            for kr = 1:nr
                restraints.depth(depth_poi).site{kr} = control.directives(d).block{kr,1};
                restraints.depth(depth_poi).r(kr) = str2double(control.directives(d).block{kr,2});
                restraints.depth(depth_poi).sigr(kr) = str2double(control.directives(d).block{kr,3});
                for karg = 1:args
                    fprintf(logfid,'  %s',control.directives(d).block{kr,karg});
                end
                fprintf(logfid,'\n');
            end
            fprintf(logfid,'\n\n');
        otherwise
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_flex:unknown_directive',...
                'directive %s is unknown',lower(control.directives(d).name));
            record_exception(exceptions{warnings},logfid);

    end
end

restraints.ddr = restraints.ddr(1:ddr_poi);
restraints.oligomer = restraints.oligomer(1:oligomer_poi);
restraints.depth = restraints.depth(1:depth_poi);

% determine residue range and initialize restrain variable
if ~isfield(restraints,'initial') || ~isfield(restraints,'final')...
        || ~isfield(restraints,'sequence')
    warnings = warnings + 1;
    exceptions{warnings} = MException('module_flex:no_sequence_directive',...
        'sequence or residue range not specified');
    record_exception(exceptions{warning},logfid);
    failed = true;
    return
end

sequence = restraints.sequence;

nres = restraints.final - restraints.initial + 1;
restrain(nres).secondary = 0;
restrain(nres).aprop = 0;
restrain(nres).bprop = 0;
restrain(nres).cprop = 0;
restrain(nres).label = [];
restrain(nres).r_beacon = [];
restrain(nres).r_intern = [];
restrain(nres).oligomer = [];
restrain(nres).depth = [];

% store propensities
% prolong propensity vectors if required
if restraints.final > maxlen
    aprop = zeros(1,restraints.final);
    aprop(1:maxlen) = restraints.aprop;
    restraints.aprop = aprop;
    bprop = zeros(1,restraints.final);
    bprop(1:maxlen) = restraints.bprop;
    restraints.bprop = bprop;
    cprop = zeros(1,restraints.final);
    cprop(1:maxlen) = restraints.cprop;
    restraints.cprop = cprop;
end
% transfer propensity information
for kres = restraints.initial:restraints.final
    krel = kres - restraints.initial + 1; % index in chain segment
    restrain(krel).secondary = 0; % legacy information 'hardcoded' secondary structure
    restrain(krel).aprop = restraints.aprop(kres);
    restrain(krel).bprop = restraints.bprop(kres);
    restrain(krel).cprop = restraints.cprop(kres);
end

% get entity coordinates
if ~isempty(entity)
    entity_xyz = entity.xyz;
else
    entity_xyz = [];
end
% get anchor information
anchorC = [];
anchorCn = [];
anchorN = [];
anchorNp = [];
anchor_chain = 'A';
% C-terminal anchor
if isfield(restraints,'c_anchor') && ~isempty(restraints.c_anchor)
    % determine anchor chain, if specified
    poi1 = strfind(restraints.c_anchor,'(');
    poi2 = strfind(restraints.c_anchor,')');
    if ~isempty(poi1) && ~isempty(poi2)
        anchor_chain = restraints.c_anchor(poi1+1:poi2-1);
    end
    anchorC = zeros(4,3);
    failed = false;
    % extract backbone coordinates of C-terminal anchor residue
    argout = get_atom(entity,'xyz',[restraints.c_anchor '.N']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorC(1,:) = argout{1};
    end
    argout = get_atom(entity,'xyz',[restraints.c_anchor '.CA']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorC(2,:) = argout{1};
    end
    argout = get_atom(entity,'xyz',[restraints.c_anchor '.C']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorC(3,:) = argout{1};
    end
    argout = get_atom(entity,'xyz',[restraints.c_anchor '.O']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorC(4,:) = argout{1};
    end
    % extract backbone coordinates of the residue after the C-trminal
    % anchor
    poi = strfind(restraints.c_anchor,')');
    if ~isempty(poi)
        resstr = restraints.c_anchor(poi+1:end);
        prefix = restraints.c_anchor(1:poi);
    else
        resstr = restraints.c_anchor;
        prefix = '';
    end
    next_res = str2double(resstr) + 1;
    c_anchor_next = sprintf('%s%i',prefix,next_res);
    anchorCn = zeros(4,3);
    argout = get_atom(entity,'xyz',[c_anchor_next '.N']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorCn(1,:) = argout{1};
    end
    argout = get_atom(entity,'xyz',[c_anchor_next '.CA']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorCn(2,:) = argout{1};
    end
    argout = get_atom(entity,'xyz',[c_anchor_next '.C']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorCn(3,:) = argout{1};
    end
    argout = get_atom(entity,'xyz',[c_anchor_next '.O']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorCn(4,:) = argout{1};
    end
    % add C anchor residue to sequence
    argout = get_residue(entity,'info',restraints.c_anchor); 
    slc = tlc2slc(argout{1}.tlc);
    sequence = [sequence slc];
    c_anchored = true;
else
    c_anchored = false;
end

if failed
    warnings = warnings + 1;
    exceptions{warnings} = MException('module_flex:c_anchor_unresolved',...
        'not all backbone atoms found for C-terminal anchor %s (aborted)',...
        restraints.c_anchor);
    record_exception(exceptions{warning},logfid);
    return
end

% N-terminal anchor
if isfield(restraints,'n_anchor') && ~isempty(restraints.n_anchor)
    % determine anchor chain, if specified
    poi1 = strfind(restraints.n_anchor,'(');
    poi2 = strfind(restraints.n_anchor,')');
    if ~isempty(poi1) && ~isempty(poi2)
        anchor_chain = restraints.n_anchor(poi1+1:poi2-1);
    end
    % determine backbone coordinates of N-terminal anchor residue
    anchorN = zeros(4,3);
    failed = false;
    argout = get_atom(entity,'xyz',[restraints.n_anchor '.N']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorN(1,:) = argout{1};
    end
    argout = get_atom(entity,'xyz',[restraints.n_anchor '.CA']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorN(2,:) = argout{1};
    end
    argout = get_atom(entity,'xyz',[restraints.n_anchor '.C']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorN(3,:) = argout{1};
    end
    argout = get_atom(entity,'xyz',[restraints.n_anchor '.O']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorN(4,:) = argout{1};
    end
    poi = strfind(restraints.n_anchor,')');
    if ~isempty(poi)
        resstr = restraints.n_anchor(poi+1:end);
        prefix = restraints.n_anchor(1:poi);
    else
        resstr = restraints.n_anchor;
        prefix = '';
    end
    % determine backbone coordinates of residue before N-terminal anchor
    % resiude
    previous_res = str2double(resstr) - 1;
    n_anchor_previous = sprintf('%s%i',prefix,previous_res);
    anchorNp = zeros(4,3);
    argout = get_atom(entity,'xyz',[n_anchor_previous '.N']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorNp(1,:) = argout{1};
    end
    argout = get_atom(entity,'xyz',[n_anchor_previous '.CA']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorNp(2,:) = argout{1};
    end
    argout = get_atom(entity,'xyz',[n_anchor_previous '.C']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorNp(3,:) = argout{1};
    end
    argout = get_atom(entity,'xyz',[n_anchor_previous '.O']);
    if isempty(argout) || isempty(argout{1})
        failed = true;
    else
        anchorNp(4,:) = argout{1};
    end
    % add N anchor residue to sequence
    argout = get_residue(entity,'info',restraints.n_anchor); 
    slc = tlc2slc(argout{1}.tlc);
    sequence = [slc sequence];
    n_anchored = true;
else
    n_anchored = false;
end

if failed
    warnings = warnings + 1;
    exceptions{warnings} = MException('module_flex:n_anchor_unresolved',...
        'not all backbone atoms found for N-terminal anchor %s (aborted)',...
        restraints.n_anchor);
    record_exception(exceptions{warning},logfid);
    return
end

n_restraints = 0;

% process distance distribution restraints
for ddr_poi = 1:length(restraints.ddr)
    label1 = restraints.ddr(ddr_poi).labels{1};
    label2 = restraints.ddr(ddr_poi).labels{2};
    curr_M_max = restraints.ddr(ddr_poi).M_max;
    curr_f_update = restraints.ddr(ddr_poi).f_update;
    for kr = 1:length(restraints.ddr(ddr_poi).site1)
        % establish type
        % ### handling of full distance distributions needs to be
        % considered later ###
        r = restraints.ddr(ddr_poi).r(kr);
        sigr = restraints.ddr(ddr_poi).sigr(kr);
        if isempty(r) || isempty(sigr)
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_flex:ddr_empty',...
                'missing distance or standard deviation for ddr %s/%s (ignored)',...
                restraints.ddr(ddr_poi).site1{kr},restraints.ddr(ddr_poi).site2{kr});
            record_exception(exceptions{warning},logfid);
            continue
        end
        if r == 0 && sigr == 0
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_flex:ddr_ignored',...
                'zero distance and standard deviation for ddr %s/%s (ignored)',...
                restraints.ddr(ddr_poi).site1{kr},restraints.ddr(ddr_poi).site2{kr});
            record_exception(exceptions{warning},logfid);
            continue
        end
        n_restraints = n_restraints+1;
        if ~isempty(restraints.ddr(ddr_poi).file{kr})
            restraint_type = 'rejection';
        elseif sigr < 0 % interpreted as lower/upper bound
            restraint_type = 'bounds';
            r = abs(r);
            sigr = abs(sigr);
        else
            restraint_type = 'Gaussian';
        end
        % try to compute label coordinates at both sites
        [argsout,entity] = get_label(entity,label1,'positions',restraints.ddr(ddr_poi).site1{kr});
        if ~isempty(argsout) && ~isempty(argsout{1})
            pos1 = argsout{1};
            [argsout,entity] = get_label(entity,label1,'populations',restraints.ddr(ddr_poi).site1{kr});
            pop1 = argsout{1};
            kres1 = 0;
        else % determine residue index in chain segment
            site1 = restraints.ddr(ddr_poi).site1{kr};
            poi = strfind(site1,')');
            if ~isempty(poi)
                site1 = site1(poi+1:end);
            end
            kres1 = str2double(site1) - restraints.initial + 1;
        end
        [argsout,entity] = get_label(entity,label2,'positions',restraints.ddr(ddr_poi).site2{kr});        
        if ~isempty(argsout) && ~isempty(argsout{1})
            pos2 = argsout{1};
            [argsout,entity] = get_label(entity,label2,'populations',restraints.ddr(ddr_poi).site2{kr});
            pop2 = argsout{1};
            kres2 = 0;
        else % determine residue index in chain segment
            site2 = restraints.ddr(ddr_poi).site2{kr};
            poi = strfind(site2,')');
            if ~isempty(poi)
                site2 = site2(poi+1:end);
            end
            kres2 = str2double(site2) - restraints.initial + 1;
        end
        % load sampling and restraint distribution, if provided
        if ~isempty(restraints.ddr(ddr_poi).file(kr)) 
            frname = restraints.ddr(ddr_poi).file(kr);
            frname = frname{1};
            if ~isempty(frname)
                data = load(frname);
                r_axis = data(:,1);
                samples = data(:,2);
                target = data(:,3);
                M = max((target-0.01*max(target))./(samples+0.01*max(samples)));
                if M > curr_M_max
                    M = curr_M_max;
                end
            else
                r_axis = [];
                samples = [];
                target = [];
                M = [];
            end
        end
        % determine type of restraint
        if kres1 > 0 && kres1 <= nres
            internal1 = true;
            mean_pos = get_relative_label(label1);
            restrain(kres1).label = mean_pos;
        else
            internal1 = false;
        end
        if kres2 > 0 && kres2 <= nres
            internal2 = true;
            mean_pos = get_relative_label(label2);
            restrain(kres2).label = mean_pos;
        else
            internal2 = false;
        end
        if ~internal1 && ~internal2 % restraint outside modelled segment is ignored
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_flex:ddr_outside_section',...
                'for ddr, neither %s nor %s is inside modelled chain segment (ignored)',...
                restraints.ddr(ddr_poi).site1{kres1},restraints.ddr(ddr_poi).site2{kres2});
            record_exception(exceptions{warning},logfid);
            continue
        end
        if internal1 && internal2 % internal segment
            kres = kres1;
            site = kres2;
            if site > kres
                kres = kres2;
                site = kres1;
            end
            if ~isfield(restrain(kres),'r_intern') % this residue does not yet have an internal restraint
                restrain(kres).r_intern(1).site = site;
                restrain(kres).r_intern(1).type = restraint_type;
                restrain(kres).r_intern(1).par1 = r;
                restrain(kres).r_intern(1).par2 = sigr;
                restrain(kres).r_intern(1).r_axis = r_axis;
                restrain(kres).r_intern(1).samples = samples;
                restrain(kres).r_intern(1).target = target;
                restrain(kres).r_intern(1).M = M;
                restrain(kres).r_intern(1).M_max = curr_M_max;
                restrain(kres).r_intern(1).f_update = curr_f_update;
            else
                kint = length(restrain(kres).r_intern) + 1;
                restrain(kres).r_intern(kint).site = site;
                restrain(kres).r_intern(kint).type = restraint_type;
                restrain(kres).r_intern(kint).par1 = r;
                restrain(kres).r_intern(kint).par2 = sigr;
                restrain(kres).r_intern(kint).r_axis = r_axis;
                restrain(kres).r_intern(kint).samples = samples;
                restrain(kres).r_intern(kint).target = target;
                restrain(kres).r_intern(kint).M = M;
                restrain(kres).r_intern(kint).M_max = curr_M_max;
                restrain(kres).r_intern(kint).f_update = curr_f_update;
            end
            continue
        end
        if internal1 % first site is in modelled segment
            kres = kres1;
            xyz = pop2'*pos2;
        else % second site is in modelled segment
            kres = kres2;
            xyz = pop1'*pos1;
        end
        if ~isfield(restrain(kres),'r_beacon') % this residue does not yet have a beacon restraint
            restrain(kres).r_beacon(1).xyz = xyz;
            restrain(kres).r_beacon(1).type = restraint_type;
            restrain(kres).r_beacon(1).par1 = r;
            restrain(kres).r_beacon(1).par2 = sigr;
            restrain(kres).r_beacon(1).r_axis = r_axis;
            restrain(kres).r_beacon(1).samples = samples;
            restrain(kres).r_beacon(1).target = target;
            restrain(kres).r_beacon(1).M = M;
            restrain(kres).r_beacon(1).M_max = curr_M_max;
            restrain(kres).r_beacon(1).f_update = curr_f_update;
        else
            kint = length(restrain(kres).r_beacon) + 1;
            restrain(kres).r_beacon(kint).xyz = xyz;
            restrain(kres).r_beacon(kint).type = restraint_type;
            restrain(kres).r_beacon(kint).par1 = r;
            restrain(kres).r_beacon(kint).par2 = sigr;
            restrain(kres).r_beacon(kint).r_axis = r_axis;
            restrain(kres).r_beacon(kint).samples = samples;
            restrain(kres).r_beacon(kint).target = target;
            restrain(kres).r_beacon(kint).M = M;
            restrain(kres).r_beacon(kint).M_max = curr_M_max;
            restrain(kres).r_beacon(kint).f_update = curr_f_update;
        end
    end
end

% process oligomer restraints
for oli_poi = 1:length(restraints.oligomer)
    label = restraints.oligomer(oli_poi).label;
    curr_M_max = restraints.oligomer(oli_poi).M_max;
    curr_f_update = restraints.oligomer(oli_poi).f_update;
    label_rel_coor = get_relative_label(label);
    multiplicity = restraints.oligomer(oli_poi).multiplicity;
    for kr = 1:length(restraints.oligomer(oli_poi).site)
        % establish type
        % ### handling of full distance distributions needs to be
        % considered later ###
        r = restraints.oligomer(oli_poi).r(kr);
        sigr = restraints.oligomer(oli_poi).sigr(kr);
        if isempty(r) || isempty(sigr)
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_flex:oligomer_empty',...
                'missing distance or standard deviation for oligomer restraint %s (ignored)',...
                restraints.oligomers(oli_poi).site{kr});
            record_exception(exceptions{warning},logfid);
            continue
        end
        if r == 0 && sigr == 0
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_flex:oligomer_ignored',...
                'zero distance and standard deviation for oligomer restraint %s (ignored)',...
                restraints.oligomer(oli_poi).site{kr});
            record_exception(exceptions{warning},logfid);
            continue
        end
        n_restraints = n_restraints+1;
        if ~isempty(restraints.ddr(ddr_poi).file{kr})
            restraint_type = 'rejection';
        elseif sigr < 0 % interpreted as lower/upper bound
            restraint_type = 'bounds';
            r = abs(r);
            sigr = abs(sigr);
        else
            restraint_type = 'Gaussian';
        end
        site = restraints.oligomer(oli_poi).site{kr};
        poi = strfind(site,')');
        if ~isempty(poi)
            site = site(poi+1:end);
        end
        kres = str2double(site) - restraints.initial + 1;
        % load sampling and restraint distribution, if provided
        if ~isempty(restraints.oligomer(oli_poi).file(kr)) 
            frname = restraints.oligomer(oli_poi).file(kr);
            frname = frname{1};
            if ~isempty(frname)
                data = load(frname);
                r_axis = data(:,1);
                samples = data(:,2);
                target = data(:,3);
                M = max((target-0.01*max(target))./(samples+0.01*max(samples)));
                if M > curr_M_max
                    M = curr_M_max;
                end
            else
                r_axis = [];
                samples = [];
                target = [];
                M = [];
            end
        end
        % check if site is inside modelled segment
        if kres > 0 && kres <= nres
            restrain(kres).label = label_rel_coor;
            if ~isfield(restrain(kres),'oligomer') % this residue does not yet have an oligomer restraint
                restrain(kres).oligomer(1).n = multiplicity;
                restrain(kres).oligomer(1).type = restraint_type;
                restrain(kres).oligomer(1).par1 = r;
                restrain(kres).oligomer(1).par2 = sigr;
                restrain(kres).oligomer(1).r_axis = r_axis;
                restrain(kres).oligomer(1).samples = samples;
                restrain(kres).oligomer(1).target = target;
                restrain(kres).oligomer(1).M = M;
                restrain(kres).oligomer(1).M_max = curr_M_max;
                restrain(kres).oligomer(1).f_update = curr_f_update;
            else
                kint = length(restrain(kres).oligomer) + 1;
                restrain(kres).oligomer(kint).n = multiplicity;
                restrain(kres).oligomer(kint).type = restraint_type;
                restrain(kres).oligomer(kint).par1 = r;
                restrain(kres).oligomer(kint).par2 = sigr;
                restrain(kres).oligomer(kint).r_axis = r_axis;
                restrain(kres).oligomer(kint).samples = samples;
                restrain(kres).oligomer(kint).target = target;
                restrain(kres).oligomer(kint).M = M;
                restrain(kres).oligomer(kint).M_max = curr_M_max;
                restrain(kres).oligomer(kint).f_update = curr_f_update;
            end 
        else
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_flex:oligomer_restraint_outside_section',...
                'site %s for oligomer restraint is not inside modelled chain segment (ignored)',...
                restraints.oligomer(oli_poi).site);
            record_exception(exceptions{warning},logfid);
            continue
        end
    end
end

% process depth restraints
for depth_poi = 1:length(restraints.depth)
    label = restraints.depth(depth_poi).label;
    if strcmpi(label,'CA') % check whether the restraint is direct to CA
        restraint_site = 'CA';
    else
        restraint_site = 'label';
        label_rel_coor = get_relative_label(label);
    end
    for kr = 1:length(restraints.depth(depth_poi).site)
        % establish type
        r = restraints.depth(depth_poi).r(kr);
        sigr = restraints.depth(depth_poi).sigr(kr);
        if isempty(r) || isempty(sigr)
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_flex:depth_empty',...
                'missing distance or standard deviation for depth restraint %s (ignored)',...
                restraints.depth(depth_poi).site{kr});
            record_exception(exceptions{warning},logfid);
            continue
        end
        if r == 0 && sigr == 0
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_flex:depth_ignored',...
                'zero distance and standard deviation for depth restraint %s (ignored)',...
                restraints.depth(depth_poi).site{kr});
            record_exception(exceptions{warning},logfid);
            continue
        end
        n_restraints = n_restraints+1;
        if sigr < 0 % interpreted as lower/upper bound
            restraint_type = 'bounds';
            r = abs(r);
            sigr = abs(sigr);
        else
            restraint_type = 'Gaussian';
        end
        site = restraints.depth(depth_poi).site{kr};
        poi = strfind(site,')');
        if ~isempty(poi)
            site = site(poi+1:end);
        end
        kres = str2double(site) - restraints.initial + 1;
        % check if site is inside modelled segment
        if kres > 0 && kres <= nres
            if strcmp(restraint_site,'label')
                restrain(kres).label = label_rel_coor;
            end
            if ~isfield(restrain(kres),'depth') % this residue does not yet have an oligomer restraint
                restrain(kres).depth(1).site = restraint_site;
                restrain(kres).depth(1).type = restraint_type;
                restrain(kres).depth(1).par1 = r;
                restrain(kres).depth(1).par2 = sigr;
            else
                kint = length(restrain(kres).depth) + 1;
                restrain(kres).depth(kint).site = restraint_site;
                restrain(kres).depth(kint).type = restraint_type;
                restrain(kres).depth(kint).par1 = r;
                restrain(kres).depth(kint).par2 = sigr;
            end 
        else
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_flex:depth_restraint_outside_section',...
                'site %s for depth restraint is not inside modelled chain segment (ignored)',...
                restraints.depth(depth_poi).site);
            record_exception(exceptions{warning},logfid);
            continue
        end
    end
end

defs = load('monomer_attributes.mat');
rescodes = zeros(1,length(sequence));
for k = 1:length(sequence)
    rescodes(k) = strfind(defs.monomers.aa_slc,sequence(k));
end

defs = load('Ramachandran_disordered.mat');

Rama_res.me = defs.Ramachandran.me;
Rama_res.ephi = defs.Ramachandran.ephi;
Rama_res.epsi = defs.Ramachandran.epsi;
Rama_res.allowed_P = defs.Ramachandran.allowed_P;
Rama_res.allowed_G = defs.Ramachandran.allowed_G;
Rama_res.allowed_gen = defs.Ramachandran.allowed_gen;

if ~isempty(control.options) || ~isempty(control.options{1})
    pmodel = str2double(control.options{1});
else
    pmodel = 0.5;
end
pthr = exp(-erfinv(pmodel)^2);
min_prob = pthr^n_restraints;

if length(control.options) >= 2
    max_models = str2double(control.options{2});
else
    max_models = 1;
end

if length(control.options) >= 3
    max_time = str2double(control.options{3});
else
    max_time = 1;
end

ntrials = 20000000; % number of Monte Carlo trials
max_seconds = 3600*max_time; % maximum runtime in seconds

free_standing = true;
reverse = false;

if n_anchored || c_anchored
    free_standing = false;
end

if c_anchored && ~n_anchored 
    reverse = true;
end

if c_anchored && n_anchored 
    closed_loop = true;
else
    closed_loop = false;
end

res1 = restraints.initial;
resend = restraints.final;
        
if ~free_standing
    prot_coor = entity_xyz;
else
    prot_coor = [];
end

success = 0;
err_count=zeros(1,11);
Ram_fixed = 0;
Ram_fix_clash = 0;
resax = res1:resend;
res_stat = zeros(1,length(resax));

fprintf(logfid,'Saving models to %s_m#.pdb\n',fname);

p_coor = cell(opt.parnum);
p_restrain = cell(opt.parnum);
p_errcode = zeros(1,opt.parnum);
p_cumprob = zeros(1,opt.parnum);
p_kres = zeros(1,opt.parnum);

run_options.min_prob = min_prob;
save_options.pdbid = pdbid;
save_options.chainIDs{1,1} = 'Z';
save_options.chainIDs{1,2} = chainid;


k_MC = 0; % number of Monte Carlo trials
tpm = NaN;

for res = 1:length(restrain)
    % checke whetehr there are beacon restraints
    if isfield(restrain(res),'r_beacon')
        % loop over all beacon restraints of this residue
        for kr = 1:length(restrain(res).r_beacon)
            % test whether we do have a full distribution restraint here
            if ~isempty(restrain(res).r_beacon(kr).samples)
                restrain(res).r_beacon(kr).updated_samples = zeros(length(restrain(res).r_beacon(kr).samples),1);
            end
        end
    end
    if isfield(restrain(res),'r_intern')
        % loop over all internal restraints of this residue
        for kr = 1:length(restrain(res).r_intern)
            % test whether we do have a full distribution restraint here
            if ~isempty(restrain(res).r_intern(kr).samples)
                restrain(res).r_intern(kr).updated_samples = zeros(length(restrain(res).r_intern(kr).samples),1);
            end
        end
    end
    if isfield(restrain(res),'oligomer')
        % loop over all oligomer restraints of this residue
        for kr = 1:length(restrain(res).oligomer)
            % test whether we do have a full distribution restraint here
            if ~isempty(restrain(res).oligomer(kr).samples)
                restrain(res).oligomer(kr).updated_samples = zeros(length(restrain(res).oligomer(kr).samples),1);
            end
        end
    end
end


tic,
while 1
    parfor kp = 1:opt.parnum % parfor
        if reverse
            [coor,errcode,restrain1,cumprob,kres] = loop_model_reverse(sequence, anchorC, anchorCn, prot_coor, restrain, Rama_res, rescodes, n_restraints, run_options);
        else
            if rejection_sampling
                [coor,errcode,restrain1,cumprob,kres] = loop_model_rs(sequence, anchorN, anchorC, anchorNp, anchorCn, prot_coor, restrain, Rama_res, rescodes, n_restraints, run_options);
            else
                [coor,errcode,restrain1,cumprob,kres] = loop_model(sequence, anchorN, anchorC, anchorNp, anchorCn, prot_coor, restrain, Rama_res, rescodes, n_restraints, run_options);
            end
            kres = kres-1;
        end
        p_coor{kp} = coor;
        p_errcode(kp) = errcode;
        p_restrain{kp} = restrain1; 
        p_cumprob(kp) = cumprob; 
        p_kres(kp) = kres;
    end
    % initialize new sampling updates
    for res = 1:length(restrain)
        % checke whetehr there are beacon restraints
        if isfield(restrain(res),'r_beacon')
            % loop over all beacon restraints of this residue
            for kr = 1:length(restrain(res).r_beacon)
                % test whether we do have a full distribution restraint here
                if ~isempty(restrain(res).r_beacon(kr).samples)
                    restrain(res).r_beacon(kr).new_samples = zeros(length(restrain(res).r_beacon(kr).samples),1);
                end
            end
        end
        if isfield(restrain(res),'r_intern')
            % loop over all internal restraints of this residue
            for kr = 1:length(restrain(res).r_intern)
                % test whether we do have a full distribution restraint here
                if ~isempty(restrain(res).r_intern(kr).samples)
                    restrain(res).r_intern(kr).new_samples = zeros(length(restrain(res).r_intern(kr).samples),1);
                end
            end
        end
        if isfield(restrain(res),'oligomer')
            % loop over all oligomer restraints of this residue
            for kr = 1:length(restrain(res).oligomer)
                % test whether we do have a full distribution restraint here
                if ~isempty(restrain(res).oligomer(kr).samples)
                    restrain(res).oligomer(kr).new_samples = zeros(length(restrain(res).oligomer(kr).samples),1);
                end
            end
        end
    end

    for kp = 1:opt.parnum
        % collect all sampling distributions
        % initialize new sampling distributions
        output_restrain = p_restrain{kp};
        for res = 1:length(restrain)
            % checke whether there are beacon restraints
            if isfield(restrain(res),'r_beacon')
                % loop over all beacon restraints of this residue
                for kr = 1:length(restrain(res).r_beacon)
                    % test whether we do have a full distribution restraint here
                    if ~isempty(restrain(res).r_beacon(kr).samples)
                        restrain(res).r_beacon(kr).updated_samples = ...
                            restrain(res).r_beacon(kr).updated_samples + output_restrain(res).r_beacon(kr).new_samples;
                        restrain(res).r_beacon(kr).new_samples = ...
                            restrain(res).r_beacon(kr).new_samples + output_restrain(res).r_beacon(kr).new_samples;
                    end
                end
            end
            if isfield(restrain(res),'r_intern')
                % loop over all internal restraints of this residue
                for kr = 1:length(restrain(res).r_intern)
                    % test whether we do have a full distribution restraint here
                    if ~isempty(restrain(res).r_intern(kr).samples)
                        restrain(res).r_intern(kr).updated_samples = ...
                            restrain(res).r_intern(kr).updated_samples + output_restrain(res).r_intern(kr).new_samples;
                        restrain(res).r_intern(kr).new_samples = ...
                            restrain(res).r_intern(kr).new_samples + output_restrain(res).r_intern(kr).new_samples;
                    end
                end
            end
            if isfield(restrain(res),'oligomer')
                % loop over all oligomer restraints of this residue
                for kr = 1:length(restrain(res).oligomer)
                    % test whether we do have a full distribution restraint here
                    if ~isempty(restrain(res).oligomer(kr).samples)
                        restrain(res).oligomer(kr).updated_samples = ...
                            restrain(res).oligomer(kr).updated_samples + output_restrain(res).oligomer(kr).new_samples;  
                        restrain(res).oligomer(kr).new_samples = ...
                            restrain(res).oligomer(kr).new_samples + output_restrain(res).oligomer(kr).new_samples;  
                    end
                end
            end
        end

        k_MC = k_MC + 1;
        coor = p_coor{kp};
        errcode = p_errcode(kp);
        kres = p_kres(kp);
        res_stat(kres) = res_stat(kres)+1;
        runtime = toc;
        if errcode == -1
            Ram_fixed = Ram_fixed + 1;
            errcode = 0;
        end
        if errcode == -4
            Ram_fixed = Ram_fixed + 1;
            Ram_fix_clash = Ram_fix_clash + 1;
            errcode = 4;
        end
        err_count(errcode+1) = err_count(errcode+1) + 1;

        if ~errcode
            tpm = runtime/err_count(1);
            success = success + 1;
            if success == 1
                bb0 = coor;
            elseif free_standing
                [rms,coor] = rmsd_superimpose(bb0,coor);
                fprintf(logfid,'Model superimposes onto first model with rmsd of %4.1f Angstroem\n',rms);
            end
             loopname = write_pdb_backbone(coor,restraints.sequence,fname,success,res1);
             pmodel = make_SCWRL4_sidegroups(loopname);
             if ~isempty(pmodel)
                 delete(loopname);
             else
                 fprintf(logfid,'Warning: Model %s could not be decorated with sidegroups\n',loopname);
                 pmodel = loopname;
             end
             [pclash,iclash] = check_decorated_loop(pmodel,prot_coor,res1,resend);

            if pclash
                err_count(10) = err_count(10) + 1;
                success = success - 1;
                fprintf(logfid,'Warning: Side groups in model %s clash with protein\n',pmodel);
                if ~keep_sidegroup_clashes
                    delete(pmodel);
                end
            elseif iclash
                err_count(11) = err_count(11) + 1;
                success = success - 1; 
                fprintf(logfid,'Warning: Side groups in model %s clash internally\n',pmodel);
                if ~keep_sidegroup_clashes
                    delete(pmodel);
                end
            else
                m_entity = get_pdb(pmodel);
                fname_valid = sprintf('%s_valid_m%i',fname,success);
                if n_anchored || c_anchored
                    entity1 = entity;
                    for resnum = res1:resend
                        resname = sprintf('R%i',resnum);
                        entity1 = add_residue(entity1,anchor_chain,resnum,m_entity.Z.(resname),m_entity.xyz);
                    end
                    put_pdb(entity1,fname_valid,save_options);
                else
                    put_pdb(m_entity,fname_valid,save_options);
                end
                delete(pmodel);
            end

            if success >= max_models
                break
            end
        end
        if opt.interactive && mod(k_MC,opt.disp_update) == 0
            ftr = (1 - k_MC/ntrials)*max_seconds;
            fti = max_seconds - runtime;
            fmo = runtime*(max_models-success)/success;
            time_left = min([fti fmo ftr]);
            hours = floor(time_left/3600);
            minutes = round((time_left-3600*floor(time_left/3600))/60);
            if minutes == 60
                hours = hours + 1;
                minutes = 0;
            end
            fprintf(1,'Time per model: %8.1f s\n',tpm);
            fprintf(1,'%i h %i min estimated run time to completion\n',hours,minutes);
            fprintf(logfid,'Successful trials: %i\n',success);
            fprintf(logfid,'Successful backbone trials: %i\n',err_count(1));
            fprintf(logfid,'Restraint violations: %5.2f%%\n',100*err_count(6)/k_MC);
            fprintf(logfid,'Internal loop clashes: %5.2f%%\n',100*(err_count(3)+err_count(8))/k_MC);
            fprintf(logfid,'Clashes with protein: %5.2f%%\n',100*(err_count(5)+err_count(7))/k_MC);
            fprintf(logfid,'Internal sidegroup clashes: %5.2f%%\n',100*err_count(11)/err_count(1));
            fprintf(logfid,'Sidegroup clashes with protein: %5.2f%%\n',100*err_count(10)/err_count(1));
            if closed_loop
                fprintf(logfid,'Loop closure failures: %5.2f%%\n',100*err_count(2)/k_MC);
                fprintf(logfid,'Ramachandran mismatches: %5.2f%%\n',100*err_count(4)/k_MC);
            end

        end
    end
    
    % update sampling information
    for res = 1:length(restrain)
        % check whether there are beacon restraints
        if isfield(restrain(res),'r_beacon')
            % loop over all beacon restraints of this residue
            for kr = 1:length(restrain(res).r_beacon)
                % test whether we do have a full distribution restraint here
                if ~isempty(restrain(res).r_beacon(kr).samples)
                    f_update = restrain(res).r_beacon(kr).f_update;
                    new_samples = restrain(res).r_beacon(kr).new_samples;
                    if sum(new_samples) > 0
                        new_samples = new_samples/sum(new_samples);
                    end
                    samples = (1-f_update)*restrain(res).r_beacon(kr).samples + f_update*new_samples;
                    samples = samples/sum(samples);
                    target = restrain(res).r_beacon(kr).target;
                    M = max((target-0.01*max(target))./(samples+0.01*max(samples)));
                    if M > restrain(res).r_beacon(kr).M_max
                        M = restrain(res).r_beacon(kr).M_max;
                    end
                    restrain(res).r_beacon(kr).samples = samples;
                    restrain(res).r_beacon(kr).M = M;
                end
            end
        end
        if isfield(restrain(res),'r_intern')
            % loop over all internal restraints of this residue
            for kr = 1:length(restrain(res).r_intern)
                % test whether we do have a full distribution restraint here
                if ~isempty(restrain(res).r_intern(kr).samples)
                    f_update = restrain(res).r_intern(kr).f_update;
                    new_samples = restrain(res).r_intern(kr).new_samples;
                    if sum(new_samples) > 0
                        new_samples = new_samples/sum(new_samples);
                    end
                    samples = (1-f_update)*restrain(res).r_intern(kr).samples + f_update*new_samples;
                    samples = samples/sum(samples);
                    target = restrain(res).r_intern(kr).target;
                    M = max((target-0.01*max(target))./(samples+0.01*max(samples)));
                    if M > restrain(res).r_intern(kr).M_max
                        M = restrain(res).r_intern(kr).M_max;
                    end
                    restrain(res).r_intern(kr).samples = samples;
                    restrain(res).r_intern(kr).M = M;
                end
            end
        end
        if isfield(restrain(res),'oligomer')
            % loop over all oligomer restraints of this residue
            for kr = 1:length(restrain(res).oligomer)
                % test whether we do have a full distribution restraint here
                if ~isempty(restrain(res).oligomer(kr).samples)
                    f_update = restrain(res).oligomer(kr).f_update;
                    new_samples = restrain(res).oligomer(kr).new_samples;
                    if sum(new_samples) > 0
                        new_samples = new_samples/sum(new_samples);
                    end
                    samples = (1-f_update)*restrain(res).oligomer(kr).samples + f_update*new_samples;
                    samples = samples/sum(samples);
                    target = restrain(res).oligomer(kr).target;
                    M = max((target-0.01*max(target))./(samples+0.01*max(samples)));
                    if M > restrain(res).oligomer(kr).M_max
                        M = restrain(res).oligomer(kr).M_max;
                    end
                    restrain(res).oligomer(kr).samples = samples;
                    restrain(res).oligomer(kr).M = M;
                end
            end
        end
    end
    
    if success >= max_models
        break
    end
    if runtime >= max_seconds
        fprintf(logfid,'Warning: Computation stopped since maximum runtime of %5.2f h was reached.\n',max_seconds/3600);
        break
    end
end
fprintf(logfid,'\nFlex module has run its course:\n');
fprintf(logfid,'\nRuntime was %i s (%8.1f s per full model):\n',round(runtime),runtime/success);
fprintf(logfid,'Successful trials: %i\n',success);
fprintf(logfid,'Successful backbone trials: %i\n',err_count(1));
fprintf(logfid,'Restraint violations: %5.2f%%\n',100*err_count(6)/k_MC);
fprintf(logfid,'Internal loop clashes: %5.2f%%\n',100*(err_count(3)+err_count(8))/k_MC);
fprintf(logfid,'Clashes with protein: %5.2f%%\n',100*(err_count(5)+err_count(7))/k_MC);
fprintf(logfid,'Internal sidegroup clashes: %5.2f%%\n',100*err_count(11)/err_count(1));
fprintf(logfid,'Sidegroup clashes with protein: %5.2f%%\n',100*err_count(10)/err_count(1));
if closed_loop
    fprintf(logfid,'Loop closure failures: %5.2f%%\n',100*err_count(2)/k_MC);
    fprintf(logfid,'Ramachandran mismatches: %5.2f%%\n',100*err_count(4)/k_MC);
end

for res = 1:length(restrain)
    % check whether there are beacon restraints
    if isfield(restrain(res),'r_beacon')
        % loop over all beacon restraints of this residue
        for kr = 1:length(restrain(res).r_beacon)
            % test whether we do have a full distribution restraint here
            if ~isempty(restrain(res).r_beacon(kr).samples)
                data = [restrain(res).r_beacon(kr).r_axis...
                    restrain(res).r_beacon(kr).samples...
                    restrain(res).r_beacon(kr).updated_samples...
                    restrain(res).r_beacon(kr).target];
                outname = sprintf('updated_samples_beacon_%i_%i.dat',res,kr);
                save(outname,'data','-ascii');
            end
        end
    end
    if isfield(restrain(res),'r_intern')
        % loop over all internal restraints of this residue
        for kr = 1:length(restrain(res).r_intern)
            % test whether we do have a full distribution restraint here
            if ~isempty(restrain(res).r_intern(kr).samples)
                data = [restrain(res).r_intern(kr).r_axis...
                    restrain(res).r_intern(kr).samples...
                    restrain(res).r_intern(kr).updated_samples...
                    restrain(res).r_intern(kr).target];
                outname = sprintf('updated_samples_intern_%i_%i.dat',res,kr);
                save(outname,'data','-ascii');
            end
        end
    end
    if isfield(restrain(res),'oligomer')
        % loop over all oligomer restraints of this residue
        for kr = 1:length(restrain(res).oligomer)
            % test whether we do have a full distribution restraint here
            if ~isempty(restrain(res).oligomer(kr).samples)
                data = [restrain(res).oligomer(kr).r_axis...
                    restrain(res).oligomer(kr).samples...
                    restrain(res).oligomer(kr).updated_samples...
                    restrain(res).oligomer(kr).target];
                outname = sprintf('updated_samples_oligomer_%i_%i.dat',res,kr);
                save(outname,'data','-ascii');
            end
        end
    end
end



function mean_pos = get_relative_label(label)

mean_pos = [];
tlc = label_by_synonym(label);
fname = sprintf('rotamer_library_%s.mat',tlc);

if ~exist(fname,'file')
    return
else
    lib = load(sprintf('rotamer_library_%s.mat',tlc));
end

position_indices = lib.rot_lib.position(:,1);
position_weights = lib.rot_lib.position(:,2);
% renormalize position weights
position_weights = position_weights/sum(position_weights);
pop = lib.rot_lib.populations;
pop = pop/sum(pop);
% positions for all rotamers
pos = zeros(length(pop),3);
for r = 1:length(pop)
    coor = lib.rot_lib.rotamers(r).coor; % r-th rotamer
    position = coor(position_indices,:);
    position = position_weights'*position; % mean position coordinates for the r-th rotamer
    pos(r,:) = position; 
end
mean_pos = pop'*pos;


function record_exception(exception,logfile)

fprintf(logfile,'### flex exception %s ###\n',exception.message);

function loopname = write_pdb_backbone(coor,sequence,fname0,model,res1)

residues='ALAARGASNASPCYSGLNGLUGLYHISILELEULYSMETPHEPROSERTHRTRPTYRVAL';
oneletter='ARNDCQEGHILKMFPSTWYV';

fname = sprintf('%s_m%i',fname0,model);

loopname = [fname '.pdb'];
wfile=fopen(loopname,'w');
for k = 1:length(sequence)
    respoi=strfind(oneletter,sequence(k));
    residue=residues(1+3*(respoi-1):3+3*(respoi-1));
    N = coor(4*k-3,:);
    CA = coor(4*k-2,:);
    C = coor(4*k-1,:);
    O = coor(4*k,:);
    fprintf(wfile,'%s%5i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           N\n','ATOM  ',4*k-3,'  N   ',residue,k+res1-1,N(1),N(2),N(3));
    fprintf(wfile,'%s%5i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           C\n','ATOM  ',4*k-2,'  CA  ',residue,k+res1-1,CA(1),CA(2),CA(3));
    fprintf(wfile,'%s%5i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           C\n','ATOM  ',4*k-1,'  C   ',residue,k+res1-1,C(1),C(2),C(3));
    fprintf(wfile,'%s%5i%s%s%6i%12.3f%8.3f%8.3f  1.00  0.00           O\n','ATOM  ',4*k,'  O   ',residue,k+res1-1,O(1),O(2),O(3));
end
fclose(wfile);


function [outname,status,result] = make_SCWRL4_sidegroups(inname)
%
% Attaches or corrects sidegroups using SCWRL4
% the program SCWRL4 needs to be on the current Matlab path
%
% inname   name of the PDB file to which sidegroups should be attached
% outname  name of the output PDB file with SCWRL4 sidegroups
%

poi = strfind(inname,'.pdb');
outname = [inname(1:poi-1) '_SCWRL4.pdb'];

s = which('scwrl4.exe');
if isempty(s)
    outname = '';
    return
end
cmd = [s ' -i ' inname ' -o ' outname];
[status,result] = dos(cmd);

function [pclash,iclash,approach_prot,approach_loop] = check_decorated_loop(loopname,prot_coor,res1,resend,min_approach)

if ~exist('min_approach','var')
    min_approach = 1.2; 
end

approach_prot = -1;
approach_loop = -1;
pclash = 1;
iclash = 1;
loop_coor = zeros(5000,3);
l_res_assign = zeros(5000,3);
fid=fopen(loopname);
if fid==-1
    return;
end
poi = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if length(tline) >= 6
        record=tline(1:6);
        resnum = str2double(tline(23:26));
        if strcmpi(record,'ATOM  ') || strcmpi(record,'HETATM')
            if length(tline) < 78 || tline(78)~='H'
                if ~strcmpi(tline(14),'H') && resnum ~= res1 && resnum ~= resend
                    poi = poi +1;
                    l_res_assign(poi) = resnum;
                    loop_coor(poi,:) = [str2double(tline(31:38)) str2double(tline(39:46)) str2double(tline(47:54))];
                end
            end
        end
    end
end
fclose(fid);

loop_coor = loop_coor(1:poi,:);
l_res_assign = l_res_assign(1:poi);

pclash = 0;
iclash = 0;

[m1,~] = size(loop_coor); % get sizes of the coordinates arrays
[m2,~] = size(prot_coor);

if m2 > 0
    a2 = repmat(sum(loop_coor.^2,2),1,m2);
    b2 = repmat(sum(prot_coor.^2,2),1,m1).';
    pair_dist = sqrt(abs(a2 + b2 - 2*loop_coor*prot_coor.'));
    min_dist = min(min(pair_dist));

    approach_prot = min_dist;
    if min_dist < min_approach
       pclash = 1;
       return
    end
end

min_dist = 1e6;
% test for minimum distance within loop
for k1 = 1:poi-1
    for k2 = k1+1:poi
        if abs(l_res_assign(k1)-l_res_assign(k2))>1
            approach = norm(loop_coor(k1,:) - loop_coor(k2,:));
            if approach < min_dist
                min_dist = approach;
            end
        end
    end
end
approach_loop = min_dist;
if min_dist < min_approach
    iclash = 1;
end

