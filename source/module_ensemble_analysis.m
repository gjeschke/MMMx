function [entity,exceptions,failed] = module_ensemble_analysis(control,logfid)
%
% MODULE_ENSEMBLE_ANALYSIS    Analyses an ensemble
%
%   [entity,exceptions] = MODULE_ENSEMBLE_ANALYSIS(control,logfid)
%   Given a conformer (entity), coordinate transformations or sidechain
%   modifications are applied
%
% INPUT
% control       control structure with fields
%               .name           'ensembleanalysis', unused
%               .options        options struct
%               .directives     array of struct with directives
%               .entity_output  number of the entity to be output
% logfid        file handle for log file, defaults to console output
%
% OUTPUT
% entity        output entity
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
% failed        flag indicating whether the module failed
%
% SUPPORTED DIRECTIVES (an ordered list is processed)
%
% addpdb        add conformers by reading pdb files
% compare       compare two ensembles
% figures       format for saving figures, 'off' switches figure saving
%               off, this is the default
% flexibility   residue-specific flexibility by Ramachandran angle
%               distribution
% flory         residue-specific flexibility by Flory's characteristic
%               ratio
% getens        get ensemble file
% gettraj       get ensemble from trajectory, requires mdtraj (actually
%               mdconvert.exe on the Matlab path)
% order         order ensemble
% paircorr      pair correlation matrix
% segments      distribution of segment-wise root-mean square end-to-end
%               distances and deviation from random-coil fit
% subsample     reduce ensemble (from a trajectory) by subsampling
% superimpose   superimpose conformers
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% initialize empty output
exceptions = {[]};
entity = [];
warnings = 0;
failed = false;


% set defaults

save_figures = false; % default is to not save figures
figure_format = 'pdf';

% default output to Matlab console if no log file identifiere was provided 
if isempty(logfid)
    logfid = 1;
end

commands = cell(1,1000); % command list
ensembles = cell(1,100); % ensemble list

cmd_poi = 0; % command pointer
ensemble_poi = 0; % entity pointer
% reorganize command line arguments
for d = 1:length(control.directives)
    clear cmd
    cmd.name = lower(control.directives(d).name);
    switch lower(control.directives(d).name)
        case {'addpdb','getens'}
            cmd_poi = cmd_poi + 1;
            ensemble_poi = ensemble_poi + 1;
            cmd.input = control.directives(d).options{1};
            cmd.ensemble = ensemble_poi;
            if length(control.directives(d).options) > 1 % the entity has an explicit internal name
                ensemble_descriptor.name = control.directives(d).options{2};
            else % the internal name is derived from the order of loading entities
                ensemble_descriptor.name = sprintf('E%i',ensemble_poi);
            end
            ensembles{ensemble_poi} = ensemble_descriptor;
            commands{cmd_poi} = cmd;
        case 'figures'
            cmd_poi = cmd_poi + 1;
            cmd.extension = control.directives(d).options{1};
            commands{cmd_poi} = cmd;
        case 'flexibility'
            cmd_poi = cmd_poi + 1;
            cmd.fname = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % a selected entity is analyzed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; % flexibility analysis is performed for current entity
            end
            commands{cmd_poi} = cmd;
        case 'flory'
            cmd_poi = cmd_poi + 1;
            cmd.fname = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % a selected entity is analyzed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; % Flory analysis is performed for current entity
            end
            commands{cmd_poi} = cmd;
        case 'subsample'
            cmd_poi = cmd_poi + 1;
            cmd.stride = str2double(control.directives(d).options{1});
            if length(control.directives(d).options) > 1 % the entity has an explicit internal name
                cmd.entity = control.directives(d).options{2};
            else % use current entity
                cmd.entity = '.';
            end
            if length(control.directives(d).options) > 2 % the entity has an explicit internal name
                cmd.output = control.directives(d).options{3};
            else % use current entity
                cmd.output = '.';
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
        otherwise
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_ensembleanalysis:unknown_directive',...
                'directive %s is unknown',lower(control.directives(d).name));
            record_exception(exceptions{warnings},logfid);
    end
end

ensembles = ensembles(1:ensemble_poi);
commands = commands(1:cmd_poi);

fprintf(logfid,'\n%i commands will be executed on %i ensembles\n\n',cmd_poi,ensemble_poi);

% run the command list

for c = 1:cmd_poi
    cmd = commands{c};
    switch cmd.name
        case 'addpdb'
            ensemble_poi = cmd.ensemble;
            ensemble_descriptor = ensembles{ensemble_poi};
            ensemble_name = ensemble_descriptor.name;
            [~,~,ext] = fileparts(cmd.input);
            if isempty(ext)
                cmd.input = strcat(cmd.input,'.pdb');
            end
            added_files = dir(cmd.input); % find all files that match the pattern
            [entity,exceptions] = get_pdb(added_files(1).name);
            if ~isempty(exceptions) && ~isempty(exceptions{1})
                warnings = warnings +1;
                exceptions{warnings} = MException('module_ensembleanalysis:file_does_not_exist',...
                    'PDB file %s could not be opened',added_files(1).name);
                record_exception(exceptions{warnings},logfid);
                return
            end
            if length(added_files) > 1
                entity.populations = ones(1,length(added_files))/length(added_files);
                for cf = 2:length(added_files)
                    [entity,exceptions] = get_pdb(added_files(cf).name,[],entity);
                    if ~isempty(exceptions) && ~isempty(exceptions{1})
                        warnings = warnings +1;
                        exceptions{warnings} = MException('module_ensembleanalysis:file_does_not_exist',...
                            'PDB file %s could not be opened',added_files(cf).name);
                        record_exception(exceptions{warnings},logfid);
                        return
                    end
                end
            end
            ensemble_descriptor.entity = entity;
            fprintf(logfid,'\nCurrent ensemble is: %s\n',ensemble_name);
            ensembles = store_ensemble(ensemble_name,entity,ensembles);
        case 'getens'
            ensemble_poi = cmd.ensemble;
            ensemble_descriptor = ensembles{ensemble_poi};
            ensemble_name = ensemble_descriptor.name;
            entity = get_ensemble(cmd.input);
            fprintf(logfid,'\nCurrent ensemble is: %s\n',ensemble_name);
            ensembles = store_ensemble(ensemble_name,entity,ensembles);
        case 'figures'
            if strcmpi(cmd.extension,'off')
                save_figures = false;
            else
                save_figures = true;
                figure_format = cmd.extension;
            end
        case 'flexibility'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'flexibility cannot be analysed for entity %s, since entity is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            if save_figures
                flexibility_analysis(c_entity,cmd.fname,figure_format);
            else
                flexibility_analysis(c_entity,cmd.fname);
            end
        case 'flory'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'Flory analysis cannot be performed for entity %s, since entity is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            if save_figures
                flory_analysis(c_entity,cmd.fname,figure_format);
            else
                flory_analysis(c_entity,cmd.fname);
            end
        case 'subsample'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'tried to subsample entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            C = length(c_entity.populations);
            selection = '{1';
            conf = 1;
            while conf + cmd.stride <= C
                conf = conf + cmd.stride;
                selection = sprintf('%s,%i',selection,conf);
            end
            selection = sprintf('%s}(*)',selection);
            c_entity = select(c_entity,selection,true);
            clear save_options
            save_options.selected = true;
            save_options.pop = true;
            put_pdb(c_entity,'temporary.pdb',save_options);
            c_entity = get_pdb('temporary.pdb');
            delete('temporary.pdb');
            if strcmp(cmd.output,'.')
                entity = c_entity;
            end
            ensembles = store_ensemble(cmd.output,c_entity,ensembles);
    end
end

function entity = retrieve_ensemble(name,ensembles,logfid)

ensemble = [];
e = 0;
while isempty(ensemble) && e < length(ensembles)
    e = e + 1;
    ensemble_descriptor = ensembles{e};
    if strcmp(name,ensemble_descriptor.name)
        entity = ensemble_descriptor.entity;
    end
end
if isempty(ensemble)
    fprintf(logfid,'Ensemble %s has not been loaded.\n',name);
end

function ensembles = store_ensemble(name,entity,ensembles)

found = false;
e = 0;
while ~found && e < length(ensembles)
    e = e + 1;
    ensemble_descriptor = ensembles{e};
    if strcmp(name,ensemble_descriptor.name)
        found = true;
        ensemble_descriptor.entity = entity;
        ensembles{e} = ensemble_descriptor;
    end
end

function record_exception(exception,logfid)

fprintf(logfid,'### ensembleanalysis exception: %s ###\n',exception.message);

function [ctag,rtag] = split_chain_residue(address)

p1 = strfind(address,'(');
p2 = strfind(address,')');
ctag = address(p1+1:p2-1);
rtag = strcat('R',address(p2+1:end));

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

function entity1 = get_backbone(logfid,entity1,conformer1,address1,entity2,conformer2,address2,mode,selection)

peptide_backbone = {'N','CA','C','O'};
nucleotide_backbone = {'P','O5_','C5_','C4_','C3_','O3_'};

[chains1,residues1,atoms1,conformers1] = split_address(address1);
[chains2,residues2,atoms2,conformers2] = split_address(address2);

entity1 = []; 
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
