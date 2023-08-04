function [entity,exceptions] = entity_from_filelist(filenames,name)
%
% ENTITY_FROM_FILELIST Generate entity from a list of PDB files
%
%   [entity,exceptions] = ENTITY_FROM_FILELIST(filenames)
%   Returns an (entity) structure in MMMx:atomic representation,
%   populations are initialized uniformly
%
%   [entity,exceptions] = GET_ENSEMBLE(filenames,name)
%   Returns an (entity) structure in MMMx:atomic representation and names
%   it
%
% INPUT
% filenames cell string of file names
% name      optional name for the entity, defaults to PDB identifier of
%           first conformer if header line exists; MMMx otherwise
%
% OUTPUT
% entity       entity structure in MMMx:atomic representation
% exceptions   error message if something went wrong, entity is empty for
%              errors but not for warnings, for warnings, only the last one
%              is reported, cell containing an empty array, if no exception
%
% Remark: The file list can also contain a single file with multiple
%         conformers (models). However, if there are several input files,
%         each one must contain only a single conformer

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty outputs
exceptions{1} = [];

if length(filenames) == 1 % single PDB, possibly with MMMx population information
    entity = get_pdb(filenames{1});
else % multiple PDB files, assume uniform populations
    cnum = length(filenames);
    [entity,exceptions] = get_pdb(filenames{1});
    if ~isempty(exceptions) && ~isempty(exceptions{1})
        exceptions{1} = MException('entity_from_filelist:file_does_not_exist',...
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
            exceptions{warnings} = MException('entity_from_filelist:file_does_not_exist',...
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
    entity.elements = repmat(entity.elements,length(filenames),1);
end

if exist('name','var')
    entity.name = name;
elseif isempty(entity.name)
    entity.name = 'MMMx';
end