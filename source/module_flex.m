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

opt.parnum = 100; % number of trials performed in the parfor loop
opt.disp_update = 200; % cycles between display updates in interactive mode
maxlen = 1000; % maximum expected loop length for memory pre-allocation

restraints.ddr(1).labels{1} = '';
restraints.oligomer(1).label = '';
restraints.depth(1).label = '';
restraints.aprop = zeros(1,maxlen);
restraints.bprop = zeros(1,maxlen);
restraints.cprop = zeros(1,maxlen);

ddr_poi = 0;
oligomer_poi = 0;
depth_poi = 0;
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
        case 'parallel'
            opt.parallel = str2double(control.directives(d).options{1});
        case 'save'
            outname = control.directives(d).options{1};
        case 'c_anchor'
            restraints.c_anchor = control.directives(d).options{1};
        case 'n_anchor'
            restraints.n_anchor = control.directives(d).options{1};
        case 'a_prop'
            [nr,~] = size(control.directives(d).block);
            if nr > 0 % at least one restraint in this block
                for karg = 1:nr
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
                for karg = 1:nr
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
                for karg = 1:nr
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
                else
                    restraints.ddr(ddr_poi).r(kr) = str2double(arg3);
                    restraints.ddr(ddr_poi).sigr(kr) = str2double(control.directives(d).block{kr,4});
                end
                if args > 4
                    arg5 = control.directives(d).block{kr,5};
                    if arg5(1) == '@'
                        restraints.ddr(ddr_poi).file{kr} = arg5(2:end);
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
                else
                    restraints.oligomer(oligomer_poi).r(kr) = str2double(arg2);
                    restraints.oligomer(oligomer_poi).sigr(kr) = str2double(control.directives(d).block{kr,3});
                end
                if args > 3
                    arg4 = control.directives(d).block{kr,4};
                    if arg4(1) == '@'
                        restraints.oligomer(oligomer_poi).file{kr} = arg4(2:end);
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
% C-terminal anchor
if isfield(restraints,'c_anchor') && ~isempty(restraints.c_anchor)
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
        if sigr < 0 % interpreted as lower/upper bound
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
                restraints.ddr(ddr_poi).site1{kr},restraints.ddr(ddr_poi).site2{kr});
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
            else
                kint = length(restrain(kres).r_intern) + 1;
                restrain(kres).r_intern(kint).site = site;
                restrain(kres).r_intern(kint).type = restraint_type;
                restrain(kres).r_intern(kint).par1 = r;
                restrain(kres).r_intern(kint).par2 = sigr;
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
        else
            kint = length(restrain(kres).r_beacon) + 1;
            restrain(kres).r_beacon(kint).xyz = xyz;
            restrain(kres).r_beacon(kint).type = restraint_type;
            restrain(kres).r_beacon(kint).par1 = r;
            restrain(kres).r_beacon(kint).par2 = sigr;
        end
    end
end

% process oligomer restraints
for oli_poi = 1:length(restraints.oligomer)
    label = restraints.oligomer(oli_poi).label;
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
        if sigr < 0 % interpreted as lower/upper bound
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
        % check if site is inside modelled segment
        if kres > 0 && kres <= nres
            restrain(kres).label = label_rel_coor;
            if ~isfield(restrain(kres),'oligomer') % this residue does not yet have an oligomer restraint
                restrain(kres).oligomer(1).n = multiplicity;
                restrain(kres).oligomer(1).type = restraint_type;
                restrain(kres).oligomer(1).par1 = r;
                restrain(kres).oligomer(1).par2 = sigr;
            else
                kint = length(restrain(kres).oligomer) + 1;
                restrain(kres).oligomer(kint).n = multiplicity;
                restrain(kres).oligomer(kint).type = restraint_type;
                restrain(kres).oligomer(kint).par1 = r;
                restrain(kres).oligomer(kint).par2 = sigr;
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

pmodel = str2double(control.options{1});
pthr = exp(-erfinv(pmodel)^2);
min_prob = pthr^n_restraints;


save MMMx_hnRNPA1_restrain restrain control sequence Rama_res anchorC anchorCn anchorN anchorNp entity_xyz rescodes min_prob n_restraints


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