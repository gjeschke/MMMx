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
% entity        input entity, can be missing if a getpdb directive exists
%
% OUTPUT
% entity        output entity
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
% failed        flag indicating whether the module failed
%
% SUPPORTED DIRECTIVES (an ordered list is processed)
%
% bilayer       determine optimal lipid bilayer plane and thickness,
% center        center a structure, also with respect to part of it
% chains        restrict an entity to only selected chains
% conformer     restrict an entity to one of its conformers
% deselenate    replaces selenocysteine and selenomethionine by their
%               equivalents cysteine and methionine  
% getpdb        load PDB file, several files can be loaded and named
% merge         build a combined structure from parts of other structures
% remove        remove residues from an entity
% renumber      renumber residues in a chain
% repack        repack all sidechains using SCWRL4
% repair        repair incomplete sidechains by replacing them with SCWRL4
% replace       replace a chain in one structure by a chain of another
%               structure
% save          save one of the loaded entities (by default the current
%               entity)
% superimpose   superimpose one structure onto another one
% symmetry      put structure into (pseudo)-symmetry frame
%               (optional) transformation of structure into bilayer frame
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
% reorganize command line arguments
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
            entity_descriptor.conformer = conformer;
            entities{entity_poi} = entity_descriptor;
            commands{cmd_poi} = cmd;
        case 'getcyana'
            cmd_poi = cmd_poi + 1;
            entity_poi = entity_poi + 1;
            cmd.input = control.directives(d).options{1};
            cmd.entity = entity_poi;
            options.maxch = 26;
            if length(control.directives(d).options) > 1 % the entity has an explicit internal name
                entity_descriptor.name = control.directives(d).options{2};
            else % the internal name is derived from the order of loading entities
                entity_descriptor.name = sprintf('E%i',entity_poi);
            end
            if length(control.directives(d).options) > 2 % a maximum number of chains is processed
                options.maxch = str2double(control.directives(d).options{3});
            end
            entity_descriptor.entity = get_cyana_pdb(cmd.input,options); % load the entity
            entity_descriptor.conformer = conformer;
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
        case 'deselenate'
            cmd_poi = cmd_poi + 1;
            if ~isempty(control.directives(d).options) && ~isempty(control.directives(d).options{1}) % a selected entity is deselenated
                cmd.entity = control.directives(d).options{1};
            else
                cmd.entity = '.'; % the current entity is centered
            end
            commands{cmd_poi} = cmd;
        case 'repair'
            cmd_poi = cmd_poi + 1;
            if ~isempty(control.directives(d).options) && ~isempty(control.directives(d).options{1}) % a selected entity is deselenated
                cmd.entity = control.directives(d).options{1};
            else
                cmd.entity = '.'; % the current entity is repaired
            end
            commands{cmd_poi} = cmd;
        case 'repack'
            cmd_poi = cmd_poi + 1;
            if ~isempty(control.directives(d).options) && ~isempty(control.directives(d).options{1}) % a selected entity is deselenated
                cmd.entity = control.directives(d).options{1};
            else
                cmd.entity = '.'; % the current entity is repacked
            end
            commands{cmd_poi} = cmd;
        case 'renumber'
            cmd_poi = cmd_poi + 1;
            cmd.address = strip(control.directives(d).options{1}); 
            cmd.shift = str2double(control.directives(d).options{2}); 
            if length(control.directives(d).options) > 2 % a chain in a selected entity is renumbered
                cmd.entity = control.directives(d).options{3};
            else
                cmd.entity = '.'; % a chain in the current entity is renumbered
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
        case 'merge'
            cmd_poi = cmd_poi + 1;
            cmd.entity = control.directives(d).options{1};
            [n,~] = size(control.directives(d).block);
            cmd.parts(n).entity = '';
            cmd.parts(n).conformer = [];
            cmd.parts(n).chain = '';
            cmd.parts(n).range = [];
            % store selections
            for k = 1:n
                cmd.parts(k).entity = control.directives(d).block{k,1};
                address = control.directives(d).block{k,2};
                [conf,chain,range] = split_conf_chain_range(address);
                cmd.parts(k).conformer = conf;
                cmd.parts(k).chain = chain;
                cmd.parts(k).range = range;
            end
            commands{cmd_poi} = cmd;        
        case 'mutate'
            cmd_poi = cmd_poi + 1;
            if ~isempty(control.directives(d).options) && ~isempty(control.directives(d).options{1})
                cmd.entity = control.directives(d).options{1};
            else
                cmd.entity = '.'; % mutation is performed in current entity
            end
            [n,~] = size(control.directives(d).block);
            cmd.mutations = cell(1,n);
            % store list of mutations
            for k = 1:n
                address = control.directives(d).block{k,1};
                [ctag,rtag] = split_chain_residue(address);
                new_aa = control.directives(d).block{k,2};
                if length(new_aa) == 3
                    slc = tlc2slc(new_aa);
                elseif length(new_aa) == 1
                    slc = upper(new_aa);
                else
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('module_prepare:unknown_residue',...
                        'mutation to %s is ignored (not a valid amino acid identifier)',new_aa);
                    record_exception(exceptions{warnings},logfid); 
                    slc = '';
                end
                cmd.mutations{k} = sprintf('%s.%s.%s',ctag,rtag,slc);
            end
            commands{cmd_poi} = cmd;        
        case 'bilayer'
            cmd_poi = cmd_poi + 1;
            cmd.type = control.directives(d).options{1};
            cmd.transform = true; % by default, coordinates are transformed
            cmd.oriented = false; % by default, bilayer orientation is also fitted
            if ~isempty(control.directives(d).options)...
                    && ~isempty(control.directives(d).options{2})...
                    && strcmpi(control.directives(d).options{2},'oriented')
                cmd.oriented = true;
            end
            if length(control.directives(d).options) > 2 % bilayer computation for a selected entity
                cmd.entity = control.directives(d).options{3};
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
            if length(control.directives(d).options) > 3 % only a selection is saved
                cmd.selection = control.directives(d).options{4};
            else
                cmd.selection = '*'; % this means that everything is saved
            end
            commands{cmd_poi} = cmd;
        case 'chains'
            cmd_poi = cmd_poi + 1;
            all_chains = char(zeros(1,length(control.directives(d).options{1})));
            chains = control.directives(d).options{1};
            cpoi = strfind(chains,')');
            chains = chains(1:cpoi-1);
            cpoi = 0;
            for k = 1:length(chains)
                if isstrprop(chains(k),'upper')
                    cpoi = cpoi + 1;
                    all_chains(cpoi) = chains(k);
                end
            end
            cmd.chains = all_chains(1:cpoi);
            if length(control.directives(d).options) > 1
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.';
            end
            commands{cmd_poi} = cmd;
        case 'conformer'
            cmd_poi = cmd_poi + 1;
            cmd.conformer = round(str2double(control.directives(d).options{1})); 
            if length(control.directives(d).options) > 1
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.';
            end
            commands{cmd_poi} = cmd;
        case 'remove'
            cmd_poi = cmd_poi + 1;
            cmd.remove = control.directives(d).options{1};
            if length(control.directives(d).options) > 1
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.';
            end
            commands{cmd_poi} = cmd;
        case 'superimpose'
            cmd_poi = cmd_poi + 1;
            cmd.mode = 'backbone'; % by default, backbone atoms are superimposed
            cmd.selection = 'direct'; % by default, no alignment is performed
            split_argument = split(control.directives(d).options{1},'.');
            cmd.moving = split_argument{1};
            cmd.moving_address = split_argument{2};
            split_argument = split(control.directives(d).options{2},'.');
            cmd.template = split_argument{1};
            cmd.template_address = split_argument{2};
            for k = 3:length(control.directives(d).options)
                switch control.directives(d).options{k}
                    case 'align'
                        cmd.selection = 'align';
                    case 'backbone'
                        cmd.mode = 'backbone';
                    case 'CA'
                        cmd.mode = 'CA';
                    case 'C4'''
                        cmd.mode = 'C4_';
                    case 'all'
                        cmd.mode = 'all';
                end
            end
            commands{cmd_poi} = cmd;
        case 'replace'
            cmd_poi = cmd_poi + 1;
            split_argument = split(control.directives(d).options{1},'.');
            cmd.original = split_argument{1};
            cmd.original_chain = strip(split_argument{2});
            split_argument = split(control.directives(d).options{2},'.');
            cmd.substitute = split_argument{1};
            cmd.substitute_chain = strip(split_argument{2});
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

entity_name = 'input';

% run the command list

for c = 1:cmd_poi
    cmd = commands{c};
    switch cmd.name
        case {'getpdb','getcyana'}
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
        case 'deselenate'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = get_entity(cmd.entity,entities,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to deselenate entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
                entity_name = cmd.entity;
            end
            [c_entity,num_Se] = deselenate(c_entity);
            fprintf(logfid,'\n%i selenium atoms were replaced by sulfur and bond lengths were corrected\n',num_Se);
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            entities = put_entity(cmd.entity,c_entity,entities);
        case 'symmetry'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                [c_entity,conformer] = get_entity(cmd.entity,entities,logfid);
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
        case 'bilayer'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                [c_entity,conformer] = get_entity(cmd.entity,entities,logfid);
                if isempty(c_entity)
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to compute bilayer for entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
                entity_name = cmd.entity;
            end
            [entity1,error] = bilayer_model(c_entity,conformer,cmd.type,cmd.oriented,logfid);
            if ~isempty(error)
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_prepare:bilayer_not_computed',...
                    'bilayer could bot be computed: %s',error);
                record_exception(exceptions{warnings},logfid);
                return
            end
            if cmd.transform
                c_entity = entity1;
                if strcmp(cmd.entity,'.')
                    entity = c_entity;
                end
                entities = put_entity(cmd.entity,c_entity,entities);
            end
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
            clear save_options
            if cmd.selection ~= '*'
                save_options.selected = true;
                c_entity = select(c_entity,cmd.selection,true);
            end
            if ~isempty(cmd.pdbid)
                save_options.pdbid = cmd.pdbid;
            end
            if exist('save_options','var')
                put_pdb(c_entity,cmd.output,save_options);
            else
                put_pdb(c_entity,cmd.output);
            end
        case 'chains'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = get_entity(cmd.entity,entities,logfid);
                if isempty(c_entity)
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to set conformer of entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            for k = 1:length(cmd.chains)
                address = sprintf('(%s)',cmd.chains(k));
                if k == 1
                    % overwrite old selection in entity, if any
                    [c_entity,exceptions] = select(c_entity,address,true);
                else
                    [c_entity,exceptions] = select(c_entity,address);
                end
            end
            save_options.selected = true;
            save_options.bfactor = true;
            save_options.pdbid = c_entity.name(1:4);
            put_pdb(c_entity,'temp.pdb',save_options);
            % reload template PDB, a bit inefficient, but safe
            c_entity = get_pdb('temp.pdb');
            delete('temp.pdb');
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            entities = put_entity(cmd.entity,c_entity,entities);            
        case 'conformer'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = get_entity(cmd.entity,entities,logfid);
                if isempty(c_entity)
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to set conformer of entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            c_entity.conformer = cmd.conformer;
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            entities = put_entity(cmd.entity,c_entity,entities);
        case 'replace'
            if strcmp(cmd.original,'.')
                c_entity = entity;
            else
                c_entity = get_entity(cmd.original,entities,logfid);
                if isempty(c_entity)
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to do substitution in entity %s, which is unknown',cmd.original);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            if strcmp(cmd.substitute,'.')
                s_entity = entity;
            else
                s_entity = get_entity(cmd.substitute,entities,logfid);
                if isempty(s_entity)
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to extract chain for substitution from entity %s, which is unknown',cmd.substitute);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            c_entity.(cmd.original_chain) = s_entity.(cmd.substitute_chain);
            [c_nat,~] = size(c_entity.xyz);
            c_entity.xyz = [c_entity.xyz; s_entity.xyz];
            c_entity.elements = [c_entity.elements, s_entity.elements];
            c_entity.occupancies = [c_entity.occupancies; s_entity.occupancies];
            c_entity = update_atom_indices(c_entity,cmd.original_chain,s_entity.(cmd.substitute_chain),c_nat);
            save_options.selected = false;
            save_options.bfactor = true;
            save_options.pdbid = c_entity.name(1:4);
            put_pdb(c_entity,'temp.pdb',save_options);
            % reload template PDB, a bit inefficient, but safe
            c_entity = get_pdb('temp.pdb');
            delete('temp.pdb');
            if strcmp(cmd.original,'.')
                entity = c_entity;
            end
            entities = put_entity(cmd.original,c_entity,entities);            
        case 'remove'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = get_entity(cmd.entity,entities,logfid);
                if isempty(c_entity)
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to remove residue %s from entity %s, which is unknown',cmd.remove,cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            [ctag,rtag] = split_chain_residue(cmd.remove);
            c_entity.(ctag) = rmfield(c_entity.(ctag),rtag);
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            entities = put_entity(cmd.entity,c_entity,entities); 
        case 'renumber'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = get_entity(cmd.entity,entities,logfid);
                if isempty(c_entity)
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to renumber residues in entity %s, which is unknown',cmd.remove,cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            chain = strip(cmd.address);
            residues = fieldnames(c_entity.(chain));
            for kr = 1:length(residues)
                residue = residues{kr};
                if strcmp(residue(1),'R') % these are residue fields
                    new_number = str2double(residue(2:end)) + cmd.shift;
                    new_residue = sprintf('R%i',new_number);
                    new_chain.(new_residue) = c_entity.(chain).(residue); 
                else
                    new_chain.(residue) = c_entity.(chain).(residue);
                end
            end
            c_entity = rmfield(c_entity,chain);
            c_entity.(chain) = new_chain; 
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            entities = put_entity(cmd.entity,c_entity,entities);     
        case 'merge'
            fprintf(logfid,'\nGenerating new entity %s by merging parts\n',cmd.entity);
            [c_entity,failed] = merge_parts(entities,cmd.parts,logfid);
            if ~failed
                entity_descriptor.name = cmd.entity; 
                entity_descriptor.entity = c_entity; 
                entity_descriptor.conformer = 1;
                entity_poi = entity_poi + 1;
                entities{entity_poi} = entity_descriptor;
            else
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_prepare:entity_unknown',...
                    'tried to merge from an entity, which is unknown');
                record_exception(exceptions{warnings},logfid);
                return                
            end
        case 'superimpose'
            if strcmp(cmd.moving,'.')
                c_entity = entity;
                c_conformer = conformer;
            else
                [c_entity,c_conformer] = get_entity(cmd.moving,entities,logfid);
                if isempty(c_entity)
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to move entity %s, which is unknown',cmd.moving);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            if strcmp(cmd.template,'.')
                t_entity = entity;
                t_conformer = conformer;
            else
                [t_entity,t_conformer] = get_entity(cmd.template,entities,logfid);
                if isempty(t_entity)
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to use entity %s as superposition template, which is unknown',cmd.template);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            c_entity = superimpose_entity(logfid,c_entity,c_conformer,cmd.moving_address,t_entity,t_conformer,cmd.template_address,cmd.mode,cmd.selection);
            if strcmp(cmd.moving,'.')
                entity = c_entity;
            end
            entities = put_entity(cmd.moving,c_entity,entities);  
        case 'repair'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = get_entity(cmd.entity,entities,logfid);
                if isempty(c_entity)
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to repair sidechains in entity %s, which is unknown',cmd.remove,cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            [c_entity,repaired] = modify_sidechains(c_entity);
            fprintf(logfid,'%i sidechains were repaired in entity %s\n',repaired,cmd.entity);
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            entities = put_entity(cmd.entity,c_entity,entities);     
        case 'repack'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = get_entity(cmd.entity,entities,logfid);
                if isempty(c_entity)
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to repack sidechains in entity %s, which is unknown',cmd.remove,cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            [c_entity,repacked] = modify_sidechains(c_entity,true);
            fprintf(logfid,'%i sidechains were repacked in entity %s\n',repacked,cmd.entity);
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            entities = put_entity(cmd.entity,c_entity,entities);     
        case 'mutate'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = get_entity(cmd.entity,entities,logfid);
                if isempty(c_entity)
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('module_prepare:entity_unknown',...
                        'tried to mutate in entity %s, which is unknown',cmd.remove,cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            [c_entity,mutated] = modify_sidechains(c_entity,false,cmd.mutations);
            fprintf(logfid,'%i sidechains were mutated in entity %s\n',mutated,cmd.entity);
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            entities = put_entity(cmd.entity,c_entity,entities);     
    end
end

function [entity,conformer] = get_entity(name,entities,logfid)

entity = [];
e = 0;
while isempty(entity) && e < length(entities)
    e = e + 1;
    entity_descriptor = entities{e};
    if strcmp(name,entity_descriptor.name)
        entity = entity_descriptor.entity;
        conformer = entity_descriptor.conformer;
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
        entities{e} = entity_descriptor;
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

function [ctag,rtag] = split_chain_residue(address)

p1 = strfind(address,'(');
p2 = strfind(address,')');
ctag = address(p1+1:p2-1);
rtag = strcat('R',address(p2+1:end));

function entity = update_atom_indices(entity,chain,s_chain,offset)

residues = fieldnames(entity.(chain));
for kr = 1:length(residues) % expand over all residues
    residue = residues{kr};
    if strcmp(residue(1),'R') % these are residue fields
        atoms = fieldnames(entity.(chain).(residue));
        for ka = 1:length(atoms) % expand over all atoms
            atom = atoms{ka};
            if isstrprop(atom(1),'upper') % these are atom fields
                indices = s_chain.(residue).(atom).tab_indices + offset;
                entity.(chain).(residue).(atom).tab_indices = indices;
            end
        end
    end
end

function chain = strip(chain_address)

pa = strfind(chain_address,'(');
if isempty(pa)
    pa = 0;
end
pe = strfind(chain_address,')');
if isempty(pe)
    pe = length(chain_address)+1;
end
chain = chain_address(pa+1:pe-1);

function entity1 = superimpose_entity(logfid,entity1,conformer1,address1,entity2,conformer2,address2,mode,selection)

peptide_backbone = {'N','CA','C','O'};
nucleotide_backbone = {'P','O5_','C5_','C4_','C3_','O3_'};

[chains1,residues1,atoms1,conformers1] = split_address(address1);
[chains2,residues2,atoms2,conformers2] = split_address(address2);

% if no conformers are addressed, use default conformers
if isempty(conformers1)
    conformers1 = conformer1;
end
if isempty(conformers2)
    conformers2 = conformer2;
end

switch mode
    case 'CA'
        atoms1 = {'CA'};
        atoms2 = atoms1;
    case 'C4_'
        atoms1 = {'C4_'};
        atoms2 = atoms1;
end

if length(chains1) ~= length(chains2) % number of addressed chains must match, otherwise correspondence cannot be established
    fprintf(logfid,'\nSuperposition failed because %i chains were addressed in template and %i chains in the moving entity\n',length(chains2),length(chains1));
    return
end

% allow for superimposing all members of an ensemble to a single conformer
% in the template
if length(conformers1) > 1 && length(conformers2) == 1
    conformers2 = conformers2*ones(size(conformers1));
end

if length(conformers1) ~= length(conformers2) % number of addressed conformers must match
    fprintf(logfid,'\nSuperposition failed because %i conformers were addressed in template and %i conformers in the moving entity\n',length(conformers2),length(conformers1));
    return
end

for conf = 1:length(conformers1)
    
    correspondence = zeros(50000,3);
    corr_poi = 0;
    % loop over all chain pairs
    for c = 1:length(chains1)
        chain1 = chains1(c);
        chain2 = chains2(c);
        % if mode is align, we now need to align the two sequences
        if strcmpi(selection(1:5),'align')
            residue_correspondence = align_in_entity(entity1,chain1,entity2,chain2);
        else
            if length(residues1) ~= length(residues2) % without alignment, both residue ranges must have the same length
                return
            else
                residue_correspondence = [residues1' residues2'];
            end
        end
        residues = fieldnames(entity1.(chain1));
        for kr = 1:length(residues) % expand over all residues
            residue1 = residues{kr};
            if strcmp(residue1(1),'R') % these are residue fields
                if strcmpi(mode,'backbone') % for backbone mode, get correct atom lists
                    slc = tlc2slc(entity1.(chain1).(residue1).name); % check whether residue is an amino acid
                    if ~isempty(slc)
                        atoms1 = peptide_backbone;
                        atoms2 = peptide_backbone;
                    else % otherwise assume that it is a nucleic acid
                        atoms1 = nucleotide_backbone;
                        atoms2 = nucleotide_backbone;
                    end
                end
                resnum1 = str2double(residue1(2:end));
                resnum2 = residue_correspondence(residue_correspondence(:,1) == resnum1,2);
                residue2 = sprintf('R%i',resnum2);
                if isfield(entity2.(chain2),residue2) % check whether residue exists
                    % check whether the two residues have the same type
                    match = strcmpi(entity1.(chain1).(residue1).name,entity2.(chain2).(residue2).name);
                    % correspondence is signalled if the residue types match or strict alignment is not requested
                    if match || ~strcmpi(selection,'align!')
                        atoms = fieldnames(entity1.(chain1).(residue1));
                        for ka = 1:length(atoms) % expand over all atoms
                            atom = atoms{ka};
                            if isstrprop(atom(1),'upper') % these are atom fields
                                % check if this is an adressed atom
                                match_atom = any(strcmp(atoms1,atom)) & any(strcmp(atoms2,atom));
                                % if atm is selected or there is no atom
                                % selection at all
                                if match_atom || (isempty(atoms1) && isempty(atoms2))
                                    if isfield(entity2.(chain2).(residue2),atom) % check that it indeed exists in second entity
                                        corr_poi = corr_poi + 1;
                                        correspondence(corr_poi,1) = entity1.(chain1).(residue1).(atom).tab_indices(conformers1(conf));
                                        correspondence(corr_poi,2) = entity2.(chain2).(residue2).(atom).tab_indices(conformers2(conf));
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    correspondence = correspondence(1:corr_poi,:);
    
    if corr_poi == 0
        fprintf(logfid,'\nSuperposition failed because no matching atoms were found\n');
        return
    end
    if corr_poi < 3
        return
    end
    coor1 = entity1.xyz(correspondence(:,1),:); % this is the moving entity
    coor2 = entity2.xyz(correspondence(:,2),:); % this is the template entity
    if corr_poi == 3
        [rmsd,~,transmat] = superimpose_3points(coor2,coor1);
    else
        [rmsd,~,transmat] = rmsd_superimpose(coor2,coor1);
    end
    xyz = entity1.xyz(entity1.index_array(:,4) == conformers1(conf),:);
    xyz = affine_coor_set(xyz,transmat);
    entity1.xyz(entity1.index_array(:,4) == conformers1(conf),:) = xyz;
    fprintf(logfid,'\nSuperposition of %i atoms in conformer %i with rmsd of %4.1f %s succeeded\n',corr_poi,conformers1(conf),rmsd,char(197));
end

function [conf,chain,range] = split_conf_chain_range(address)

range = [];

conf = 1;
pa = strfind(address,'{');
pe = strfind(address,'}');
if ~isempty(pa) && ~isempty(pe)
    conf = str2double(address(pa+1:pe-1));
    if pa > 1
        pre = address(1:pa-1);
    else
        pre = '';
    end
    if pe < length(address)
        past = address(pe+1:end);
    else
        past = '';
    end
    address = [pre past];
end

poia = strfind(address,'(');
poie = strfind(address,')');
if isempty(poie)
    chain = '';
else
    chain = address(poia+1:poie-1);
    address = address(poie+1:end);
end
if strcmp(chain,'*')
    range = [1,1e6];
    return
end
residues = split(address,'-');
if ~isempty(residues{1})
    range(1) = str2double(residues{1});
end
if ~isempty(residues{2})
    range(2) = str2double(residues{2});
end
if length(range) == 1
    range(2) = range(1);
end

function [entity,failed] = merge_parts(entities,parts,logfid)

failed = false;

entity.name = 'MMMx';
entity.populations = 1;
entity.selected = 0;
entity.index_array = zeros(100000,5,'uint16');
xyz = zeros(100000,3);
elements = zeros(1,100000,'uint8');
occupancies = ones(100000,1,'uint8');
atnum = 0;
for p = 1:length(parts)
    ctag = char(double('A')+p-1); % chain tag for this part in new entity
    entity.(ctag).selected = false;
    entity.(ctag).index = p;
    conf = parts(p).conformer;
    chain = parts(p).chain;
    range = parts(p).range(1):parts(p).range(2);
    c_entity = get_entity(parts(p).entity,entities,logfid);
    if isempty(c_entity)
        fprintf(logfid,'\nERROR: Entity %s could not be retrieved.\n',parts(p).entity);
        failed = true;
        entity = [];
        return
    else
        fprintf(logfid,'New chain %s derived from [%s](%s)%i-%i\n',ctag,c_entity.name,chain,parts(p).range);
    end
    for r = range % loop over all selected residues
        residue = sprintf('R%i',r);
        entity.(ctag).(residue).index = 1 + r - range(1);
        entity.(ctag).(residue).selected = false;
        entity.(ctag).(residue).selected_rotamers = 1;
        entity.(ctag).(residue).name = c_entity.(chain).(residue).name;
        entity.(ctag).(residue).locations = ' ';
        entity.(ctag).(residue).populations = 1;
        atoms = fieldnames(c_entity.(chain).(residue));
        for ka = 1:length(atoms)
            atom = atoms{ka};
            if isstrprop(atom(1),'upper') % these are atom fields
                this_atom = c_entity.(chain).(residue).(atom);
                at_index = this_atom.tab_indices(conf);
                atnum = atnum + 1;
                this_atom.tab_indices = atnum;
                xyz(atnum,:) = c_entity.xyz(at_index,:);
                elements(atnum) = c_entity.elements(at_index);
                entity.(ctag).(residue).(atom) = this_atom;
            end
        end
    end
end
entity.xyz = xyz(1:atnum,:);
entity.elements = elements(1:atnum);
entity.occupancies = occupancies(1:atnum);

