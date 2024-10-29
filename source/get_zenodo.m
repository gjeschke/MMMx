function [entity,exceptions] = get_zenodo(ZenodoID,fname,options)
%
% GET_ZENODO Load ensemble structure from Zenodo
%
%   [entity,exceptions] = GET_ZENODO(ZenodoID,fname,options)
%   Returns an (entity) structure in MMMx:atomic representation
%
% INPUT
% ZenodoID      Zenodo dataset identifier, for instance '6384003' for 
%               the hnRNP A1 dataset, empty output and exception, when this 
%               Zenodo ID does not exist
% fname         file name in the Zenodo dataset, for instance 
%               'hnRNPA1_unrestrained_raw_ensemble.zip', empty output and
%               exception, when this file does not exist or cannot be
%               interpreted as an ensemble
% options       .keep_file  if true, the downloaded and extracted file
%                           is kept, defaults to false (file is deleted)
%
% OUTPUT
% entity       entity structure in MMMx:atomic representation 
%              .name            Z, followed by first three digits of Zenodo ID
%              in addition to PDB entities, PED entities contain
%              .origin          'Zenodo:%i.%s', where %i is the Zenodo 
%                               identifier and %s the file name 
% exceptions   error message if something went wrong, entity is empty for
%              errors but not for warnings, for warnings, only the last one
%              is reported, cell containing an empty array, if no exception
%
% REMARKS
% - The Zenodo file can be a single PDB file or a ZIP archive or a (gzipped)
%   TAR archive. 
% - For a single PDB file, models are taken as conformers.
% - If an archive contains several PDB files, they are taken as conformers
% - An archive may additionally contain an MMMx ensemble (.ens) file. In
%   this case, the ensemble is constructed according to this file
% - A gzipped input file must contain only a single tar archive

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

% initialize empty outputs
exceptions{1} = [];
entity = [];

if ~exist('options','var') || isempty(options)
    options.keep_file = false;
end

query = sprintf('https://zenodo.org/api/records/%s',ZenodoID);
try
    zenodo_info = webread(query);
catch exception
    exceptions{1} = exception;
    return;
end

query = '';
k = 0;
while k < length(zenodo_info.files)
    k = k + 1;
    if strcmpi(zenodo_info.files(k).key,fname)
        query = zenodo_info.files(k).links.self;
        break
    end
end

try
    websave(fname,query);
catch exception
    exceptions{1} = exception;
    return;
end

% resolve file format cases 
[~,~,ext] = fileparts(fname);

switch ext
    case '.ens'
        filenames{1} = fname;
    case '.pdb'
        filenames{1} = fname;
    case '.zip'
        filenames = unzip(fname);
    case '.gz'
        filenames = gunzip(fname);
        if length(filenames) == 1
            [~,~,ext2] = fileparts(filenames{1});
            if strcmpi(ext2,'.tar')
                filenames = untar(filenames{1});
            end
        end
    case '.tar'
        filenames = untar(fname);
end

% search for a (single) .ens file
ensemble_file = '';
for k = 1:length(filenames)
    [~,~,ext] = fileparts(filenames{k});
    if strcmpi(ext,'.ens')
        if isempty(ensemble_file)
            ensemble_file = filenames{k};
        else
            exceptions = {MException('get_Zenodo:more_than_one_ensemble','Zenodo package %s.%s contains more than one .ens file',ZenodoID,fname)};
            return
        end
    end
end

if ~isempty(ensemble_file)
    entity = get_ensemble(ensemble_file);
else
    poi = 0;
    for k = 1:length(filenames)
        [~,~,ext] = fileparts(filenames{k});
        if strcmpi(ext,'.pdb')
            poi = poi + 1;
            filenames{poi} = filenames{k}; %#ok<AGROW> 
        end
    end
    filenames = filenames(1:poi);
    if poi == 0
        exceptions{1} = MException('get_Zenodo:No .ens and no PDB file',...
            'Zenodo file %s does not contain ensemble information or a PDB file',fname);
        return
    end
    if length(filenames) == 1 % single PDB, possibly with MMMx population information
        entity = get_pdb(filenames{1});
    else % multiple PDB files, assume uniform populations
        cnum = length(filenames);
        [entity,exceptions] = get_pdb(filenames{1});
        if ~isempty(exceptions) && ~isempty(exceptions{1})
            exceptions{1} = MException('get_Zenodo:file_does_not_exist',...
                'PDB file %s could not be opened',filenames{1});
            return
        end
        n_atoms = length(entity.occupancies);
        occupancies = zeros(length(filenames)*n_atoms,1,'uint8');
        occupancies(1:n_atoms) = entity.occupancies;
        xyz = zeros(length(filenames)*n_atoms,3);
        xyz(1:n_atoms,:) = entity.xyz;
        index_array = zeros(length(filenames)*n_atoms,5,'uint16');
        index_array(1:n_atoms,:) = entity.index_array;
        chains = fieldnames(entity);
        for kc = 1:length(chains)
            chain = chains{kc};
            if isstrprop(chain(1),'upper') % chain fields start with a capital
                residues = fieldnames(entity.(chain));
                for kr = 1:length(residues) % expand over all residues
                    residue = residues{kr};
                    if strcmp(residue(1),'R') % these are residue fields
                        atoms = fieldnames(entity.(chain).(residue));
                        for ka = 1:length(atoms) % expand over all atoms
                            atom = atoms{ka};
                            if isstrprop(atom(1),'upper') % these are atom fields
                                at_index = entity.(chain).(residue).(atom).tab_indices(1);
                                entity.(chain).(residue).(atom).tab_indices = zeros(1,cnum);
                                entity.(chain).(residue).(atom).tab_indices(1) = at_index;
                            end
                        end
                    end
                end
            end
        end
        for cf = 2:length(filenames)
            pdb_options.conf = cf;
            pdb_options.atoff = (cf-1)*n_atoms;
            [c_entity,exceptions] = get_pdb(filenames{cf},pdb_options);
            chains = fieldnames(entity);
            for kc = 1:length(chains)
                chain = chains{kc};
                if isstrprop(chain(1),'upper') % chain fields start with a capital
                    residues = fieldnames(entity.(chain));
                    for kr = 1:length(residues) % expand over all residues
                        residue = residues{kr};
                        if strcmp(residue(1),'R') % these are residue fields
                            atoms = fieldnames(entity.(chain).(residue));
                            for ka = 1:length(atoms) % expand over all atoms
                                atom = atoms{ka};
                                if isstrprop(atom(1),'upper') % these are atom fields
                                    at_index = c_entity.(chain).(residue).(atom).tab_indices(1);
                                    entity.(chain).(residue).(atom).tab_indices(cf) = at_index;
                                end
                            end
                        end
                    end
                end
            end
            if ~isempty(exceptions) && ~isempty(exceptions{1})
                exceptions{warnings} = MException('get_Zenodo:file_does_not_exist',...
                    'PDB file %s could not be opened',filenames{cf});
                return
            end
            occupancies((cf-1)*n_atoms+1:cf*n_atoms) = c_entity.occupancies;
            xyz((cf-1)*n_atoms+1:cf*n_atoms,:) = c_entity.xyz;
            c_entity.index_array(:,4) = cf*ones(n_atoms,1,'uint16');
            index_array((cf-1)*n_atoms+1:cf*n_atoms,:) = c_entity.index_array;
        end
        entity.populations = ones(length(filenames),1)/length(filenames);
        entity.xyz = xyz;
        entity.occupancies = occupancies;
        entity.index_array = index_array;
    end
end

entity.origin = sprintf('Zenodo:%s.%s',ZenodoID,fname);
entity.name = sprintf('Z%s',ZenodoID(1:3));

if ~options.keep_file
    delete(fname);
end

