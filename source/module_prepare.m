function [entity,exceptions,failed] = module_prepare(control,logfid,entity)
%
% MODULE_PREPARE    Performs modifications of an atomic structure in order
%                   to prepare it for modelling or analysis
%
%   [entity,exceptions] = MODULE_PREPARE(control,logfid,entity)
%   Given a conformer (entity), coordinate transformations or sidechain
%   modifications are applied
%
% INPUT
% control       control structure with fields
%               .name           'prepare', unused
%               .options        options struct
%               .directives     array of struct with directives
%               .entity_output  number of the entity to be output
% logfid        file handle for log file, defaults to console output
% entity        input entity, defines rigid-body coordinates, must exist
%
% OUTPUT
% entity        output entity
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
% failed        flag indicating whether the module failed
%
% SUPPORTED DIRECTIVES (an ordered list is processed)
%
% getpdb        load PDB file, several files can be loaded and named
% center        center a structure, also with respect to part of it
% symmetry      put structure into (pseudo)-symmetry frame
% bilayer       determine optimal lipid bilayer plane and thickness,
%               (optional) transformation of structure into bilayer frame
% save          save one of the loaded entities (by default the current
%               entity)
% chains        restrict an entity to only selected chains
% conformer     restrict an entity to one of its conformers
% remove        remove residues from an entity
% sidechains    [requires SCWRL4] replace non-native side chains by native
%               equivalents, or add, or repack, or repair side chains
% superimpose   superimpose one structure onto another one
% replace       replace a chain in one structure by a chain of another
%               structure
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% initialize empty output
exceptions = {[]};
warnings = 0;
failed = false;


% set defaults

save_name = 'MMMx_prepare'; % default name for saving output entity

conformer = 1;

min_sym_axis = 5; % minimum length of a symmetry axis in Angstroem

% default output to Matlab console if no log file identifiere was provided 
if isempty(logfid)
    logfid = 1;
end

commands = cell(1,1000); % command list
entities = cell(1,100); % entity list

cmd_poi = 0; % command pointer
entity_poi = 0; % entity pointer
% read restraints
for d = 1:length(control.directives)
    clear cmd
    cmd.name = lower(control.directives(d).name);
    switch lower(control.directives(d).name)
        case 'getpdb'
            cmd_poi = cmd_poi + 1;
            entity_poi = entity_poi + 1;
            cmd.input = control.directives(d).options{1};
            cmd.entity = entity_poi;
            if length(control.directives(d).options) > 1 % the entity has an explicit internal name
                entity_descriptor.name = control.directives(d).options{2};
            else % the internal name is derived from the order of loading entities
                entity_descriptor.name = sprintf('E%i',entity_poi);
            end
            entity_descriptor.entity = get_pdb(cmd.input); % load the entity
            entities{entity_poi} = entity_descriptor;
            commands{cmd_poi} = cmd;
        case 'center'
            cmd_poi = cmd_poi + 1;
            cmd.address = '*'; % default is the whole structure
            if ~isempty(control.directives(d).options) && ~isempty(control.directives(d).options{1})
                cmd.address = control.directives(d).options{1}; % center w.r.t. part of the structure
            end
            if length(control.directives(d).options) > 1 % a selected entity is centered
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; % the current entity is centered
            end
            commands{cmd_poi} = cmd;
        case 'symmetry'
            cmd_poi = cmd_poi + 1;
            if ~isempty(control.directives(d).options) && ~isempty(control.directives(d).options{1})
                cmd.mode = control.directives(d).options{1}; % superposition mode [backbone|CA|C4'|all]
            else
                cmd.mode = 'backbone'; % by default, the backbones are superimposed 
            end
            if length(control.directives(d).options) > 1 % a selected entity transformed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; % symmetry transform is performed for current entity
            end
            [n,~] = size(control.directives(d).block);
            cmd.Cn = n; % order of symmetry axis
            cmd.parts = cell(1,n);
            % store selections in the n symmetry-equivalent parts
            for k = 1:n
                cmd.parts{k} = control.directives(d).block{k,1};
            end
            commands{cmd_poi} = cmd;
        case 'bilayer'
            cmd_poi = cmd_poi + 1;
            cmd.transform = false; % by default, coordinates are not transformed
            if ~isempty(control.directives(d).options)...
                    && ~isempty(control.directives(d).options{1})...
                    && strcmpi(control.directives(d).options{1},'transform')
                cmd.transform = true;
            end
            if length(control.directives(d).options) > 1 % bilayer computation for a selected entity
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; % bilayer is computed for current entity
            end
            commands{cmd_poi} = cmd;
        case 'save'
            cmd_poi = cmd_poi + 1;
            cmd.output = save_name; % use the default output name (can cause overwrite)
            if ~isempty(control.directives(d).options)...
                    && ~isempty(control.directives(d).options{1})
                    cmd.output = control.directives(d).options{1}; % explicit file name was assigned
            end
            if length(control.directives(d).options) > 1 % a selected entity is saved
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; % the current entity is saved
            end
            if length(control.directives(d).options) > 2 % a PDB identifier is provided
                cmd.pdbid = control.directives(d).options{3};
            else
                cmd.pdbid = ''; % this means that the original PDB identifier from the loaded structure is used
            end
            commands{cmd_poi} = cmd;
        case 'chains'
            cmd_poi = cmd_poi + 1;
            cmd.chains = char(zeros(1,length(control.directives(d).options)));
            for k = 1:length(control.directives(d).options)
                cmd.chains(k) = control.directives(d).options{k}(1);
            end
            commands{cmd_poi} = cmd;
        case 'conformer'
            cmd_poi = cmd_poi + 1;
            cmd.conformer = round(str2double(control.directives(d).options{1}));
            commands{cmd_poi} = cmd;
        case 'remove'
            cmd_poi = cmd_poi + 1;
            cmd.remove = cell(1,length(control.directives(d).options));
            for k = 1:length(control.directives(d).options)
                cmd.remove{k} = control.directives(d).options{k};
            end
            commands{cmd_poi} = cmd;
        case 'sidechains'
            cmd_poi = cmd_poi + 1;
            cmd.mode = 'repack'; % by default, all sidechains are reoptimized
            cmd.modres = false; % by default, replace applies to all non-native residues
            for k = 1:length(control.directives(d).options)
                switch control.directives(d).options{k}
                    case 'replace'
                        cmd.mode = 'replace';
                    case 'repack'
                        cmd.mode = 'repack';
                    case 'repair'
                        cmd.mode = 'repair';
                    case 'modres'
                        cmd.modres = true; % apply replace only to residues with MODRES record
                end
            end
            if strcmp(cmd.mode,'repair') % repair mode allows for selections
                [n,~] = size(control.directives(d).block);
                cmd.selected = cell(1,n);
                % extract chains, residues, maximum distances
                for k = 1:n
                    cmd.selected{k} = control.directives(d).block{k,1};
                end
            else
                cmd.selected = {};
            end
            commands{cmd_poi} = cmd;
        case 'superimpose'
            cmd_poi = cmd_poi + 1;
            cmd.mode = 'backbone'; % by default, backbone atoms are superimposed
            cmd.selection = 'align'; % by default, aligned residues are superimposed
            cmd.moving = '.';
            cmd.template = '';
            for k = 1:length(control.directives(d).options)
                split_argument = split(control.directives(d).options{k},'.');
                switch split_argument{1}
                    case 'moving'
                        cmd.moving = split_argument{2};
                    case 'template'
                        cmd.template = split_argument{2};
                    case 'align'
                        cmd.selection = 'align';
                    case 'backbone'
                        cmd.mode = 'backbone';
                    case 'CA'
                        cmd.mode = 'CA';
                    case 'C4'''
                        cmd.mode = 'C4p';
                    case 'all'
                        cmd.mode = 'all';
                    otherwise
                        cmd.selection = parameter; % superimpose only a selection
                end
            end
            commands{cmd_poi} = cmd;
        case 'replace'
            cmd_poi = cmd_poi + 1;
            split_argument = split(control.directives(d).options{1},'.');
            cmd.original = split_argument{1};
            cmd.original_chain = split_argument{2};
            split_argument = split(control.directives(d).options{2},'.');
            cmd.substitute = split_argument{1};
            cmd.substitute_chain = split_argument{2};
            commands{cmd_poi} = cmd;
        otherwise
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_prepare:unknown_directive',...
                'directive %s is unknown',lower(control.directives(d).name));
            record_exception(exceptions{warnings},logfid);
    end
end

entities = entities(1:entity_poi);
commands = commands(1:cmd_poi);

fprintf(logfid,'\n%i commands will be executed on %i structures\n\n',cmd_poi,entity_poi);

entity_poi = 0;
entity_name = 'input';

% run the command list

for c = 1:cmd_poi
    cmd = commands{c};
    switch cmd.name
        case 'getpdb'
            entity_poi = cmd.entity;
            entity_descriptor = entities{entity_poi};
            entity_name = entity_descriptor.name;
            entity = entity_descriptor.entity;
            fprintf(logfid,'Current structure is now: %s\n',entity_name);
        case 'center'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = get_entity(cmd.entity,entities,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to center entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
                entity_name = cmd.entity;
            end
            selection = 'all atoms';
            all_atoms = true;
            if ~strcmp(cmd.address,'*')
                selection = cmd.address;                
                all_atoms = false;
            end
            fprintf(logfid,'Centering structure %s with respect to %s\n',entity_name,selection);
            if all_atoms
                coor = c_entity.xyz;
            else
                coor = get_coor(c_entity,selection);
            end
            center = mean(coor,1);
            coor = c_entity.xyz;
            [n_atoms,~] = size(coor);
            c_entity.xyz = coor - repmat(center,n_atoms,1);
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            entities = put_entity(cmd.entity,c_entity,entities);
        case 'symmetry'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = get_entity(cmd.entity,entities,logfid);
                if isempty(c_entity)
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to symmetrize entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
                entity_name = cmd.entity;
            end
            [coor,error] = extract_symmetry_coordinates(c_entity,cmd.parts,cmd.mode,conformer);
            if ~isempty(error)
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_prepare:wrong_symmetry_specification',...
                    error);
                record_exception(exceptions{warnings},logfid);
                return
            end
            % compute symmetry axis
            [p0,v] = regression_line_3D(coor);
            if norm(v) < min_sym_axis
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_prepare:too_short_symmetry_axis',...
                   'Symmetry axis must have minimum length of %4.1f %s',min_sym_axis,char(197));
               record_exception(exceptions{warnings},logfid);
               return
            end
            v = v/norm(v);
            if ~isempty(p0) && ~isempty(v)
                th = acos(v(3));
                if norm(v(1:2)) > 1e-6
                    v = v(1:2)/norm(v(1:2));
                else
                    v = [0,0];
                end
                phi = atan2(v(2),v(1));
                transmat1 = affine('translation',-p0);
                transmat2 = affine('Euler',[-phi,-th,0]);
                c_entity.xyz = affine_coor_set(c_entity.xyz,{transmat1,transmat2});
            else
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_prepare:symmetry_axis_ill_defined',...
                   'Symmetry axis is ill-defined');
               record_exception(exceptions{warnings},logfid);
               return
            end
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            entities = put_entity(cmd.entity,c_entity,entities);
        case 'save'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = get_entity(cmd.entity,entities,logfid);
                if isempty(c_entity)
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to save entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
                entity_name = cmd.entity;
            end
            if isempty(cmd.pdbid)
                put_pdb(c_entity,cmd.output);
            else
                clear save_options
                save_options.pdbid = cmd.pdbid;
                put_pdb(c_entity,cmd.output,save_options);
            end
    end
end

disp('Aber hallo!');

function entity = get_entity(name,entities,logfid)

entity = [];
e = 0;
while isempty(entity) && e < length(entities)
    e = e + 1;
    entity_descriptor = entities{e};
    if strcmp(name,entity_descriptor.name)
        entity = entity_descriptor.entity;
    end
end
if isempty(entity)
    fprintf(logfid,'Entity %s has not been loaded.\n',name);
end

function entities = put_entity(name,entity,entities)

found = false;
e = 0;
while ~found && e < length(entities)
    e = e + 1;
    entity_descriptor = entities{e};
    if strcmp(name,entity_descriptor.name)
        found = true;
        entity_descriptor.entity = entity;
        entities{e} = entity_descriptor.entity;
    end
end

function [coor,error] = extract_symmetry_coordinates(entity,parts,atommode,conformer)

peptide_backbone = {'N','CA','C','O'};
nucleotide_backbone = {'P','O5_','C5_','C4_','C3_','O3_'};
    
error = '';
coor = zeros(50000,3*length(parts));
cpoi = 0;
chains = char(zeros(1,length(parts)));
res_start = -1000*ones(1,length(parts));
res_end = 1e6*ones(1,length(parts));
% extract chains and residue ranges of symmetry-equivalent parts
for p = 1:length(parts)
    address = parts{p};
    paropen = strfind(address,'(');
    parclose = strfind(address,')');
    try
        chains(p) = address(paropen+1:parclose-1); 
    catch
        error = sprintf('Wrong chain specification in %s',address);
        coor = [];
        return
    end
    if length(address) > parclose
        range = split(address(parclose+1:end),'-');
        res_start(p) = str2double(range{1});
        if length(range) < 2
            res_end(p) = res_start(p);
        else
            res_end(p) = str2double(range{2});
        end
    end
    if p == 1
        residues = 1 + res_end(1) - res_start(1);
        if isnan(residues) || residues < 1
            error = sprintf('Wrong residue range specification: %s',address);
            coor = [];
            return
        end
    else
        if residues ~= 1 + res_end(p) - res_start(p)
            error = sprintf('Inconsistent residue ranges between %s and %s',parts{1},address);
            coor = [];
            return
        end
    end
end
% extract all atom coordinates of symmetry-equivalent parts that exist in
% all parts
residues = fieldnames(entity.(chains(1)));
for kr = 1:length(residues) % expand over all residues
    residue = residues{kr};
    equivalent = cell(1,length(parts));
    if strcmp(residue(1),'R') % these are residue fields
        equivalent{1} = residue;
        resnum = str2double(residue(2:end)); % residue number
        if resnum >= res_start(1) && resnum <= res_end(1) % the residue is in the selected range
            % check whether the corresponding residue is present in all
            % parts
            all_present = true;
            for p = 2:length(parts)
                this_num = resnum + res_start(p) - res_start(1);
                restag = sprintf('R%i',this_num);
                if ~isfield(entity.(chains(p)),restag)
                    all_present = false;
                else
                    equivalent{p} = restag;
                end
            end
            if ~all_present % skip residues that do not exist in all symmetry-related parts
                continue
            end
            atoms = fieldnames(entity.(chains(1)).(residue));
            for ka = 1:length(atoms) % expand over all atoms
                atom = atoms{ka};
                if isstrprop(atom(1),'upper') % these are atom fields
                    % check whether this atom should be considered
                    switch atommode
                        case 'backbone'
                            if ~any(strcmp(peptide_backbone,atom)) && ~any(strcmp(nucleotide_backbone,atom))
                                continue
                            end
                        case 'CA'
                            if ~strcmpi(atom,'CA')
                                continue
                            end
                        case 'C4'''
                            if ~strcmpi(atom,'C4_')
                                continue
                            end
                    end
                    % check whether it exists in the other parts
                    all_present = true;
                    for p = 2:length(parts)
                        if ~isfield(entity.(chains(p)).(equivalent{p}),atom)
                            all_present = false;
                        end
                    end
                    % atom is present in all parts, record coordinates
                    if all_present
                        cpoi = cpoi + 1;
                        index = entity.(chains(1)).(restag).(atom).tab_indices(conformer);
                        coor(cpoi,1:3) = entity.xyz(index,:);
                        for p = 2:length(parts)
                            bas = 3*(p-1);
                            index = entity.(chains(p)).(equivalent{p}).(atom).tab_indices(conformer);
                            coor(cpoi,bas+1:bas+3) = entity.xyz(index,:);
                        end
                    end
                end
            end
        end
    end
end
coor0 = coor(1:cpoi,:); % restrict coordinate array to number of atoms found
if cpoi < 4
    error = 'Less than four atoms found that exist in all symmetry-equivalent parts';
    return
end
% compute center coordinates
coor = zeros(cpoi,3);
for p = 1:length(parts)
    bas = 3*(p-1);
    coor = coor + coor0(:,bas+1:bas+3);
end
coor = coor/length(parts);

function record_exception(exception,logfid)

fprintf(logfid,'### prepare exception: %s ###\n',exception.message);