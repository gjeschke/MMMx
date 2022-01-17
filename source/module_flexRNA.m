function [entity,exceptions,failed] = module_flexRNA(control,logfid,entity)
%
% MODULE_FLEXRNA Run FlexRNA module for modelling a nucleic acid loop
%
%   [exceptions,entity] = MODULE_FLEXRNA(control_file)
%   processes FlexRNA controls and restraints and runs the appropriate
%   segment modellers and stitches the RNA
%
% INPUT
% control       control structure with fields
%               .name           'flexrna', unused
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
% failed        flag indicating whether the module failed
%
% SUPPORTED DIRECTIVES
%
% sequence      nucleotide sequence, mandatory, all others are optional
% parallel      number of parallel runs per block
% expand        expand an MMMx:RigiFlex model by computing flexible linker
%               for each rigid-body arrangement
% interactive   if present, information on progress is displayed, optional
%               argument: trials between updates
% save [fn]     specify file name fn for output, defaults to mmmx_flex
% anchor_3p     3'-terminal anchor nucleotide
% anchor_5p     5'-terminal anchor nucleotide
% ddr           distance distribution restraints
% 

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% initialize empty output
exceptions = {[]};
warnings = 0;
failed = false;

initial_ensemble = '';
added_conformers = '';

if ~exist('entity','var')
    entity = [];
end
% set defaults

nent = 1; % number of entities to be processed
expand_rba = false;
use_file_list = false;

opt.parnum = 100; % number of trials performed in the parfor loop
opt.disp_update = 200; % cycles between display updates in interactive mode
opt.interactive = false;

fname = 'mmmx_flexRNA';
pdbid = 'MMMX';
chainid = 'A';
first_conformer = 1;

restraints.ddr(1).labels{1} = '';
restraints.anchor_3p = '';
restraints.anchor_5p = '';

ddr_poi = 0;
% read restraints
for d = 1:length(control.directives)
    switch lower(control.directives(d).name)
        case 'sequence'
            % first, last nucleotide number and sequence
            restraints.initial =  str2double(control.directives(d).options{1});
            restraints.final =  str2double(control.directives(d).options{2});
            restraints.sequence = control.directives(d).options{3};
        case 'interactive'
            opt.interactive = true;
            if ~isempty(control.directives(d).options) && ~isempty(control.directives(d).options{1})
                opt.disp_update = str2double(control.directives(d).options{1});
            end
        case 'expand'
            expand_rba = true;
            if ~isempty(control.directives(d).options) && ~isempty(control.directives(d).options{1})
                load(control.directives(d).options{1});
            end
        case {'initial','getpdb'}
            initial_ensemble = control.directives(d).options{1};
        case 'skipto'
            first_conformer = str2double(control.directives(d).options{1});
        case 'addpdb'
            added_conformers = control.directives(d).options{1};
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
        case 'anchor_3p'
            restraints.anchor_3p = control.directives(d).options{1};
        case 'anchor_5p'
            restraints.anchor_5p = control.directives(d).options{1};
        case 'ddr'
            ddr_poi = ddr_poi + 1; % increase ddr block counter
            restraints.ddr(ddr_poi).labels{1} = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % different labels
                restraints.ddr(ddr_poi).labels{2} = control.directives(d).options{2};
            else % same label at both sites
                restraints.ddr(ddr_poi).labels{2} = control.directives(d).options{1};
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
                exceptions{warnings} = MException('module_flexRNA:empty_ddr_block', 'ddr block %i is empty',d);
                record_exception(exceptions{warnings},logfid);
            end
            if args < 3 % misformed restraint line
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_flexRNA:misformed_ddr', 'ddr restraint has less than three arguments');
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
                    restraints.ddr(ddr_poi).file{kr} = arg3(2:end);
                    restraints.ddr(ddr_poi).r(kr) = NaN;
                    restraints.ddr(ddr_poi).sigr(kr) = NaN;
                else
                    restraints.ddr(ddr_poi).r(kr) = str2double(arg3);
                    restraints.ddr(ddr_poi).sigr(kr) = str2double(control.directives(d).block{kr,4});
                end
                if args > 4
                    arg5 = control.directives(d).block{kr,5};
                    if ~isempty(arg5) && arg5(1) == '@'
                        restraints.ddr(ddr_poi).file{kr} = arg5(2:end);
                    end
                end
            end
        otherwise
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_flexRNA:unknown_directive',...
                'directive %s is unknown',lower(control.directives(d).name));
            record_exception(exceptions{warnings},logfid);

    end
end

% determine residue range and initialize restrain variable
if ~isfield(restraints,'initial') || ~isfield(restraints,'final')...
        || ~isfield(restraints,'sequence')
    warnings = warnings + 1;
    exceptions{warnings} = MException('module_flexRNA:no_sequence_directive',...
        'sequence or residue range not specified');
    record_exception(exceptions{warning},logfid);
    failed = true;
    return
end

restraints.ddr = restraints.ddr(1:ddr_poi);

nres = restraints.final - restraints.initial + 1;
restrain(nres).label = [];
restrain(nres).r_beacon = [];
restrain(nres).r_intern = [];

fname_basis = fname;


initial_files = rd_ensemble_definition(initial_ensemble);
[~,~,ext] = fileparts(added_conformers);
if isempty(ext)
    added_conformers = strcat(added_conformers,'.pdb');
end
added_files = dir(added_conformers); % find all files that match the pattern


if expand_rba
    if isfield(entity,'rba_populations')
        nent = length(entity.rba_populations);
        entity0 = entity;
    else
        warnings = warnings + 1;
        exceptions{warnings} = MException('module_flexRNA:no_rba',...
                'Rigid-body expansion switched off, as there are no rigid bodies in input entity');   
        record_exception(exceptions{warnings},logfid);
        expand_rba = false;
    end
else
    nent = length(initial_files) + length(added_files);
    if nent > 0
        use_file_list = true;
        file_list = cell(1,nent);
        for c = 1:length(initial_files)
            file_list{c} = initial_files(c).name;
        end
        for c = 1:length(added_files)
            file_list{length(initial_files)+c} = added_files(c).name;
        end
    else
        nent = 1;
        use_file_list = false;
    end
end

for kent = first_conformer:nent
    if expand_rba
        fname = sprintf('%s_rba_%i',fname_basis,kent);
        entity = get_rba(entity0,kent);
    elseif use_file_list
        fname = sprintf('%s_conformer_%i',fname_basis,kent);
        inname = file_list{kent};
        fprintf(logfid,'Conformer %i derived from input %s\n',kent,inname);
        entity = get_pdb(inname);
    else
        fname = fname_basis;
    end
    if ~isempty(entity)
        entity = select(entity,'(*)*',true);
        if ~isempty(restraints.anchor_5p)
            entity = select(entity,restraints.anchor_5p,false,true);
        end
        if ~isempty(restraints.anchor_3p)
            entity = select(entity,restraints.anchor_3p,false,true);
        end
        environment = get_coor(entity);
        entity = select(entity,'',true);
    else
        environment = false;
    end

    anchor_3p = [];
    slc_3p = 'A';
    if ~isempty(restraints.anchor_3p)
        if isempty(entity)
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_flexRNA:anchor_needs_entity',...
                'anchored RNA loop requested, but no template provided');
            record_exception(exceptions{warnings},logfid);
        else
            [anchor_3p,slc_3p] = get_anchor_3p(entity,restraints.anchor_3p);
        end
    end
    
    anchor_5p = [];
    slc_5p = 'A';
    if ~isempty(restraints.anchor_5p)
        if isempty(entity)
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_flexRNA:anchor_needs_entity',...
                'anchored RNA loop requested, but no template provided');
            record_exception(exceptions{warnings},logfid);
        else
            [anchor_5p,slc_5p] = get_anchor_5p(entity,restraints.anchor_5p);
        end
    end

    sequence = [slc_3p restraints.sequence slc_5p];

    n_restraints = 0;
    
    % process distance distribution restraints
    for ddr_poi = 1:length(restraints.ddr)
        label1 = restraints.ddr(ddr_poi).labels{1};
        label2 = restraints.ddr(ddr_poi).labels{2};
        for kr = 1:length(restraints.ddr(ddr_poi).site1)
            % establish type
            r = restraints.ddr(ddr_poi).r(kr);
            sigr = restraints.ddr(ddr_poi).sigr(kr);
            if isempty(r) || isempty(sigr)
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_flexRNA:ddr_empty',...
                    'missing distance or standard deviation for ddr %s/%s (ignored)',...
                    restraints.ddr(ddr_poi).site1{kr},restraints.ddr(ddr_poi).site2{kr});
                record_exception(exceptions{warning},logfid);
                continue
            end
            if r == 0 && sigr == 0
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_flexRNA:ddr_ignored',...
                    'zero distance and standard deviation for ddr %s/%s (ignored)',...
                    restraints.ddr(ddr_poi).site1{kr},restraints.ddr(ddr_poi).site2{kr});
                record_exception(exceptions{warning},logfid);
                continue
            end
            n_restraints = n_restraints+1;
            if ~isempty(restraints.ddr(ddr_poi).file{kr})
                restraint_type = 'distribution';
            elseif sigr < 0 % interpreted as lower/upper bound
                restraint_type = 'bounds';
                r = abs(r);
                sigr = abs(sigr);
            else
                restraint_type = 'Gaussian';
            end
            % try to compute label coordinates at both sites, this also
            % generates the label coordinates in the entity
            [argsout,entity] = get_label(entity,label1,'positions',restraints.ddr(ddr_poi).site1{kr});
            if ~isempty(argsout) && ~isempty(argsout{1})
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
                kres2 = 0;
            else % determine residue index in chain segment
                site2 = restraints.ddr(ddr_poi).site2{kr};
                poi = strfind(site2,')');
                if ~isempty(poi)
                    site2 = site2(poi+1:end);
                end
                kres2 = str2double(site2) - restraints.initial + 1;
            end
            % load restraint distribution, if provided
            if ~isempty(restraints.ddr(ddr_poi).file(kr))
                frname = restraints.ddr(ddr_poi).file(kr);
                frname = frname{1};
                if ~isempty(frname)
                    data = load(frname);
                    r_axis = data(:,1);
                    target = data(:,2);
                    target = target/sum(target);
                    if max(r_axis) < 20
                        r_axis = 10*r_axis;
                    end
                else
                    r_axis = [];
                    target = [];
                end
            end
            % determine type of restraint
            if kres1 > 0 && kres1 <= nres
                internal1 = true;
            else
                internal1 = false;
            end
            if kres2 > 0 && kres2 <= nres
                internal2 = true;
            else
                internal2 = false;
            end
            if ~internal1 && ~internal2 % restraint outside modelled segment is ignored
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_flexRNA:ddr_outside_section',...
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
                    restrain(kres).r_intern(1).site1 = restraints.ddr(ddr_poi).site1{kr};
                    restrain(kres).r_intern(1).site2 = restraints.ddr(ddr_poi).site2{kr};
                    restrain(kres).r_intern(1).label1 = restraints.ddr(ddr_poi).labels{1};
                    restrain(kres).r_intern(1).label2 = restraints.ddr(ddr_poi).labels{2};
                    restrain(kres).r_intern(1).type = restraint_type;
                    restrain(kres).r_intern(1).par1 = r;
                    restrain(kres).r_intern(1).par2 = sigr;
                    restrain(kres).r_intern(1).r_axis = r_axis;
                    restrain(kres).r_intern(1).target = target;
                else
                    kint = length(restrain(kres).r_intern) + 1;
                    restrain(kres).r_intern(kint).site1 = restraints.ddr(ddr_poi).site1{kr};
                    restrain(kres).r_intern(kint).site2 = restraints.ddr(ddr_poi).site2{kr};
                    restrain(kres).r_intern(kint).label1 = restraints.ddr(ddr_poi).labels{1};
                    restrain(kres).r_intern(kint).label2 = restraints.ddr(ddr_poi).labels{2};
                    restrain(kres).r_intern(kint).site = site;
                    restrain(kres).r_intern(kint).type = restraint_type;
                    restrain(kres).r_intern(kint).par1 = r;
                    restrain(kres).r_intern(kint).par2 = sigr;
                    restrain(kres).r_intern(kint).r_axis = r_axis;
                    restrain(kres).r_intern(kint).target = target;
                end
                continue
            end
            if internal1 % first site is in modelled segment
                kres = kres1;
            else % second site is in modelled segment
                kres = kres2;
            end
            if ~isfield(restrain(kres),'r_beacon') % this residue does not yet have a beacon restraint
                restrain(kres).r_beacon(1).type = restraint_type;
                restrain(kres).r_beacon(1).site1 = restraints.ddr(ddr_poi).site1{kr};
                restrain(kres).r_beacon(1).site2 = restraints.ddr(ddr_poi).site2{kr};
                restrain(kres).r_beacon(1).label1 = restraints.ddr(ddr_poi).labels{1};
                restrain(kres).r_beacon(1).label2 = restraints.ddr(ddr_poi).labels{2};
                restrain(kres).r_beacon(1).par1 = r;
                restrain(kres).r_beacon(1).par2 = sigr;
                restrain(kres).r_beacon(1).r_axis = r_axis;
                restrain(kres).r_beacon(1).target = target;
            else
                kint = length(restrain(kres).r_beacon) + 1;
                restrain(kres).r_beacon(kint).type = restraint_type;
                restrain(kres).r_beacon(kint).site1 = restraints.ddr(ddr_poi).site1{kr};
                restrain(kres).r_beacon(kint).site2 = restraints.ddr(ddr_poi).site2{kr};
                restrain(kres).r_beacon(kint).label1 = restraints.ddr(ddr_poi).labels{1};
                restrain(kres).r_beacon(kint).label2 = restraints.ddr(ddr_poi).labels{2};
                restrain(kres).r_beacon(kint).par1 = r;
                restrain(kres).r_beacon(kint).par2 = sigr;
                restrain(kres).r_beacon(kint).r_axis = r_axis;
                restrain(kres).r_beacon(kint).target = target;
            end
        end
    end
    
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
    
    max_seconds = 3600*max_time; % maximum runtime in seconds
    
    if isempty(anchor_3p) && isempty(anchor_5p)
        environment = [];
    end
       
    fprintf(logfid,'Saving models to %s_m#.pdb\n',fname);
    
    save_options.pdbid = pdbid;
    
    % load the fragment library
    
    HNP = load('HNP_nucleotide_library.mat');
    acodes = [];
    transmat = [];
    ecodes = [];
    
    chain1 = '';
    chain2 = '';
    if ~isempty(restraints.anchor_5p)
        chain1  = sep_address(restraints.anchor_5p);
    end
    if ~isempty(restraints.anchor_3p)
        chain2  = sep_address(restraints.anchor_3p);
    end


    % classify the loop type, sequence always contains anchor nucleotides,
    % they are dummy nucleotides if there are no anchors
    if length(sequence) > 3 % at least two nucleotides to be modelled
        if isempty(anchor_3p) && isempty(anchor_5p)
            RNA_type = 'freestanding';
            stitch = [chainid chainid];
        elseif isempty(anchor_3p)
            RNA_type = 'trailing';
            [acodes,transmat] = get_RNA_initial_anchor(logfid,HNP.fragments,anchor_5p,upper(sequence(1)));
            if isempty(acodes)
                failed = true;
                return
            end
            stitch = [chain1 chain1];
        elseif isempty(anchor_5p)
            RNA_type = 'leading';
            ecodes = get_anchor_fragments(anchor_3p,'back',HNP.shortfrag);
            stitch = [chain2 chain2];
        else
            RNA_type = 'linker';
            [acodes,transmat] = get_RNA_initial_anchor(logfid,HNP.fragments,anchor_5p,upper(sequence(1)));
            if isempty(acodes)
                failed = true;
                return
            end
            ecodes = get_anchor_fragments(anchor_3p,'back',HNP.shortfrag);
            stitch = [chain1 chain2];
        end
    else % this is a single nucleotide
        if isempty(anchor_3p) && isempty(anchor_5p)
            RNA_type = 'any_nucleotide';
            stitch = [chainid chainid];
        elseif isempty(anchor_5p)
            RNA_type = 'leading_nucleotide';
            ecodes = get_anchor_fragments(anchor_3p,'back',HNP.shortfrag);
            stitch = [chain2 chain2];
        elseif isempty(anchor_3p)
            RNA_type = 'trailing_nucleotide';
            [acodes,transmat] = get_RNA_initial_anchor(logfid,HNP.fragments,anchor_5p,upper(sequence(1)));
            if isempty(acodes)
                failed = true;
                return
            end
            stitch = [chain1 chain1];
        else
            RNA_type = 'linking_nucleotide';
            [acodes,transmat] = get_RNA_initial_anchor(logfid,HNP.fragments,anchor_5p,upper(sequence(1)));
            if isempty(acodes)
                failed = true;
                return
            end
            ecodes = get_anchor_fragments(anchor_3p,'back',HNP.shortfrag);
            stitch = [chain1 chain2];
        end
    end
    
    ntoffset = restraints.initial - 1;
    options.maxtime = max_time;
    
    errors = zeros(3,1);
    
    success = 0;
    tstart = tic;
    err = 0;
    entity1 = [];
    
    while 1 && err ~= -3 % if anchor-anchor distance is too long, the loop finishes before maximum runtime
        for trial = 1:max_models-success
            switch RNA_type
                case {'freestanding','leading','trailing','linker','leading_nucleotide','trailing_nucleotide'}
                    [ecoor,atomtags,seq,~,~,err] = RNA_loop_model(logfid,sequence,HNP,anchor_3p,...
                        acodes,transmat,ecodes,options,ntoffset,environment);
                    if err < 0
                        errors(-err) = errors(-err) + 1;
                    end
                case 'any_nucleotide'
                    base = sequence(2);
                    code = 1 + round(rand*(length(HNP.fragments)-1));
                    coor = HNP.fragments(code).(base).coor;
                    [m,~] = size(coor);
                    atomtags = HNP.fragments(code).(base).atomtags;
                    ecoor = [(1+ntoffset)*ones(m-3,1) coor(2:m-2,:)];
                    atomtags = atomtags(2:end-2);
                    seq = sequence(2);
                    err = 0;
                case 'linking_nucleotide'
                    [ecoor,atomtags] = get_link_nt(logfid,sequence,HNP,anchor_3p,...
                        anchor_5p,ntoffset,environment);
                    seq = sequence(2);
                    err = 0;
            end
            if err == 0
                threaded = check_RNA_threading(environment,ecoor(:,2:4),atomtags);
                if isempty(threaded) || threaded
                    err = 1;
                    fprintf(logfid,'RNA is threaded through hole in protein. Continuing search.\n');
                end
            end
            if err == -3
                break
            end
        end
        if err == 0
            entity1 = insert_RNA(entity,ecoor,atomtags,seq,stitch);
            fname_temp = sprintf('%s_temp_m%i',fname,success);
            put_pdb(entity1,fname_temp,save_options);
            entity1 = get_pdb([fname_temp '.pdb']);
            prob = test_ddr(restrain,entity1,stitch);
            if prob >= min_prob
                success = success + 1;
                fname_valid = sprintf('%s_valid_m%i',fname,success);
                put_pdb(entity1,fname_valid,save_options);
                fprintf(logfid,'Model %s with restraint fulfilment probability %5.3f was saved.\n',fname_valid,prob);
            else
                fprintf(logfid,'Model rejected by too low restraint fulfilment probability (%5.3f).\n',prob);
            end
            delete([fname_temp '.pdb']);
        end
        runtime = toc(tstart);

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
            % fprintf(logfid,'Restraint violations: %5.2f%%\n',100*err_count(6)/k_MC);
        end
        
        if success >= max_models
            break
        end
        if runtime >= max_seconds
            fprintf(logfid,'Warning: Computation stopped since maximum runtime of %5.2f h was reached.\n',max_seconds/3600);
            break
        end
    end
    fprintf(logfid,'\nFlexRNA module has run its course:\n');
    fprintf(logfid,'\nRuntime was %i s (%8.1f s per full model):\n',round(runtime),runtime/success);
    fprintf(logfid,'Successful trials: %i\n\n',success);

end

entity = entity1;


function record_exception(exception,logfile)

fprintf(logfile,'### flex exception %s ###\n',exception.message);


function [anchor,slc] = get_anchor_5p(entity,address)

[anchor,slc] = get_anchor_pseudo_torsion(entity,address,true);

function [anchor,slc] = get_anchor_3p(entity,address)

[anchore,slc] = get_anchor_pseudo_torsion(entity,address);
anchor = anchore(2:4,:);
        
function [chaintag,resnum,conformer]  = sep_address(address)

conformer = [];
poi1 = strfind(address,'{');
poi2 = strfind(address,'}');
if ~isempty(poi1) && ~isempty(poi2)
    conformer = str2double(address(poi1+1:poi2-1));
    address = [address(1:poi1-1) address(poi2+1:end)];
end

poi1 = strfind(address,'(');
poi2 = strfind(address,')');
if ~isempty(poi1) && ~isempty(poi2)
    chaintag = address(poi1+1:poi2-1);
else
    chaintag = '';
    poi2 = 0;
end
resnum = str2double(address(poi2+1:end));

function [anchor,slc] = get_anchor_pseudo_torsion(entity,address,O35p)
% anchor = get_anchor_pseudo_torsion(entity,address,O35p)
% 
% retrieves a coordinate array for the pseudo-torsion defining atom
% coordinates of an anchor nucleotide, C4'(i-1), P(i), C4'(i), P(i+1),
% and, optionally, O5'(i), O3'(i)
% non-existing atoms have NaN coordinates
%
% entity    MMMx:atomic entity
% address   MMMx address of anchor nucleotide, must comprise chain tag and 
%           residue number, can contain conformer number
% O35p      optional flag that requests the O3' and O5' coordinates,
%           defaults to false
% 
% anchor    4x3 or 6x3 coordinate array
%
% G. Jeschke, 4.5.2021

[chain,resnum,conformer]  = sep_address(address);


res = sprintf('R%i',resnum);
prev_res = sprintf('R%i',resnum-1);
next_res = sprintf('R%i',resnum+1);

slc = entity.(chain).(res).name;
% treatment for DNA, remove leading 'D'
if length(slc) > 1
    slc = slc(end);
end

if exist('O35p','var') && ~isempty(O35p) && O35p
    anchor = nan(6,3);
    at_indices = entity.(chain).(res).O5_.tab_indices;
    if ~isempty(conformer)
        at_indices = at_indices(conformer);
    end
    anchor(5,:) = entity.xyz(at_indices,:);
    at_indices = entity.(chain).(res).O3_.tab_indices;
    if ~isempty(conformer)
        at_indices = at_indices(conformer);
    end
    anchor(6,:) = entity.xyz(at_indices,:);
else
    anchor = nan(4,3);
end

if isfield(entity.(chain),prev_res)
    at_indices = entity.(chain).(prev_res).C4_.tab_indices;
    if ~isempty(conformer)
        at_indices = at_indices(conformer);
    end
    anchor(1,:) = entity.xyz(at_indices,:);
end
at_indices = entity.(chain).(res).P.tab_indices;
if ~isempty(conformer)
    at_indices = at_indices(conformer);
end
anchor(2,:) = entity.xyz(at_indices,:);
at_indices = entity.(chain).(res).C4_.tab_indices;
if ~isempty(conformer)
    at_indices = at_indices(conformer);
end
anchor(3,:) = entity.xyz(at_indices,:);
if isfield(entity.(chain),next_res)
    at_indices = entity.(chain).(next_res).P.tab_indices;
    if ~isempty(conformer)
        at_indices = at_indices(conformer);
    end
    anchor(4,:) = entity.xyz(at_indices,:);
end

function entity = insert_RNA(entity,ecoor,atomtags,seq,stitch)
% adds RNA atoms to an entity
%
% entity    MMMx:atomic entity
% ecoor     extended coordinates (residue number, xyz) of RNA atoms
% atomtags  atom tags for these atoms
% seq       sequence in single-letter code
% stitch    chain identifiers for the RNA itself and for downstream RNA
%           that should be stitched to the inserted RNA
%
% G. Jeschke, 4.5.2021

if isempty(entity) || ~isfield(entity,'xyz')
    entity = [];
end

chain = stitch(1);
[nat,~] = size(ecoor);
res_offset = ecoor(1,1) - 1; % residue offset
if ~isempty(entity)
    [nat0,~] = size(entity.xyz);
    [~,cid] = get_coor(entity,sprintf('(%s)',chain));
    index_array = zeros(nat0+nat,5);
    index_array(1:nat0,:) = entity.index_array;
    cid = index_array(cid(1),1);
else
    nat0 = 0;
    index_array = zeros(nat0+nat,5);
    cid = 1;
end
xyz = zeros(nat0+nat,3);
occ = zeros(nat0+nat,1,'uint8');

if ~isempty(entity)
    xyz(1:nat0,:) = entity.xyz;
    occ0 = entity.occupancies;
    occ(1:nat0) = occ0;
else
    entity.populations = 1;
    entity.name = 'MMMx';
    entity.selected = 1;
    entity.water_selected = 0;
    entity.(stitch(1)).selected = 0; 
    entity.(stitch(1)).index = cid; 
end
resnum = -1e6;
atom_index = 0;
for kat = 1:nat
    xyz_index = nat0 + kat;
    if ecoor(kat,1) ~= resnum
        atom_index = 0;
        resnum = ecoor(kat,1);
        resstr = sprintf('R%i',resnum);
        resname = seq(ecoor(kat,1)-res_offset);
        entity.(chain).(resstr).index = ecoor(kat,1);
        entity.(chain).(resstr).selected = 0;
        entity.(chain).(resstr).selected_rotamers = 1;
        entity.(chain).(resstr).name = resname;
        entity.(chain).(resstr).locations = ' ';
        entity.(chain).(resstr).populations = 1;
    end
    atom_index = atom_index + 1;
    xyz(xyz_index,:) = ecoor(kat,2:4);
    occ(xyz_index) = 100;
    index_array(xyz_index,1) = cid;
    index_array(xyz_index,2) = ecoor(kat,1);
    index_array(xyz_index,3) = atom_index;
    index_array(xyz_index,4) = 1;
    index_array(xyz_index,5) = 1;
    atname = atomtags{kat};
    for kk = 1:length(atname)
        if atname(kk) == ''''
            atname(kk) = '_';
        end
    end
    entity.(chain).(resstr).(atname).element = atname(1);
    entity.(chain).(resstr).(atname).charge = 0;
    entity.(chain).(resstr).(atname).bfactor = 0;
    entity.(chain).(resstr).(atname).selected = 0;
    entity.(chain).(resstr).(atname).selected_locations = 1;
    entity.(chain).(resstr).(atname).tab_indices = xyz_index;
    entity.(chain).(resstr).(atname).index = atom_index;
end
entity.xyz = xyz;
entity.occupancies = occ;
entity.index_array = index_array;

% if this is a linker, rename second chain
if stitch(2) ~= stitch(1) % chain must be stitched
    residues = fieldnames(entity.(stitch(2)));
    % the following does not produce a consistent MMMx entity, but the 
    % entity can be correctly saved as PDB file
    for kres = 1:length(residues)
        resname = residues{kres};
        if resname(1) == 'R'
            entity.(stitch(1)).(resname) = entity.(stitch(2)).(resname);
        end
    end
    entity = rmfield(entity,stitch(2));
end

function prob = test_ddr(restrain,entity,stitch)

prob = 1;

nres = length(restrain);
n_ddr = 0;
for kres = 1:nres
    for kb = 1:length(restrain(kres).r_beacon)
        site1 = restrain(kres).r_beacon(kb).site1;
        site2 = restrain(kres).r_beacon(kb).site2;
        if stitch(2) ~= stitch(1)
            poi = strfind(site1,stitch(2));
            if ~isempty(poi)
                site1(poi) = stitch(1);
            end
            poi = strfind(site2,stitch(2));
            if ~isempty(poi)
                site2(poi) = stitch(1);
            end
        end
        if ~contains(site1,')')
            site1 = sprintf('(%s)%s',stitch(1),site1);
        end
        if ~contains(site2,')')
            site2 = sprintf('(%s)%s',stitch(1),site2);
        end
        [r_axis0,distribution,entity] = distance_distribution(entity,site1,...
            restrain(kres).r_beacon(kb).label1,site2,restrain(kres).r_beacon(kb).label2);
        % reject cases where the conformer cannot be labelled
        if isempty(distribution)
            prob = 0;
            return
        end
        overlap = restraint_overlap(r_axis0,distribution,restrain(kres).r_beacon(kb));        
        prob = prob*overlap;
        n_ddr = n_ddr + 1;
    end
    for kb = 1:length(restrain(kres).r_intern)
        site1 = restrain(kres).r_intern(kb).site1;
        site2 = restrain(kres).r_intern(kb).site2;
        [r_axis0,distribution,entity] = distance_distribution(entity,site1,...
            restrain(kres).r_intern(kb).label1,site2,restrain(kres).r_intern(kb).label2);
        % reject cases where the conformer cannot be labelled
        if isempty(distribution)
            prob = 0;
            return
        end
        overlap = restraint_overlap(r_axis0,distribution,restrain(kres).r_intern(kb));
        prob = prob*overlap;
        n_ddr = n_ddr + 1;
    end
end

prob = prob^(1/n_ddr);

function overlap = restraint_overlap(r_axis0,distribution,restraint)

switch restraint.type
    case 'distribution'
        r_axis = restraint.r_axis;
        target = restraint.target;
    case 'Gaussian'
        r_axis = r_axis0';
        arg = (r_axis-restraint.par1)/(sqrt(2)*restraint.par2);
        target = exp(-arg.^2);
        target = target/sum(target);
    case 'bounds'
        r_axis = r_axis0';
        target = ones(size(r_axis));
        target(r_axis<restraint.par1) = 0;
        target(r_axis>restraint.par2) = 0;
        target = target/sum(target);
 end
% figure; clf; hold on;
distribution = interp1(r_axis0,distribution,r_axis,'pchip',0);
distribution = distribution/sum(distribution);
max_distr = max(distribution);
target = max_distr*target/max(target);
% plot(r_axis,distribution);
% plot(r_axis,target);
overlap = sum(min([distribution.';target.']));