function [argout,exceptions] = get_atom(entity,attribute,address)
%
% GET_ATOM Retrieves attributes of atoms
%
%   [argout,exceptions] = GET_ATOM(entity,attribute)
%   Provides attribute values and possibly exceptions for atoms selected in
%   an entity
%
%   [argout,exceptions] = GET_ATOM(entity,attribute,address)
%   Provides attribute values and possibly exceptions for atoms selected by
%   address
%
% INPUT
% entity       entity in an MMMx format, must be provided
% address      MMMx address, 'selected' refers to the current selection
%              defaults to 'selected'
% attribute    atom attribute to be retrieved, string, defaults to
%              'info', all locations of the atom are included
%              attribute   output                           Type
%              ------------------------------------------------------------
%              bfactor     B factor                         double
%              charge      atom charge                      int
%              element     atomic number                    int8
%              info        .name          atom name         string
%                                         R# for a rotamer  # is a number
%                          .indices       index vector      (1,5) uint16
%                          .atom_index    atom array index  int
%              population  rotamer population or occupancy  double
%              xyz         Cartesian coordinates            (1,3) double
%              ecoor       extended coordinates as a        (N,5) double
%                          single array, columns
%                          1  : atomic numbers
%                          2-4: Cartesian coordinates
%                          5  : populations
%                          the second cell holds the (N,1)
%                          corresponding atom indices
%
% OUTPUT
% argout       cell array of outputs, see above
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
argout = cell(1,10000);
exceptions = {[]};

% set default arguments
if ~exist('address','var')
    address = 'selected';
end
if ~exist('attribute','var')
    attribute = 'info';
end

% select the objects by the provided address
entity = select(entity,address,true);

outputs = 0; % counter for the number of outputs
index_vector = zeros(1,5,'uint16');
conformers = entity.selected; % selected conformers
% scan entity for selected locations
chains = fieldnames(entity);
[m,~] = size(entity.index_array);
all_atom_indices = 1:m;
if strcmp(attribute,'ecoor')
    all_indices = zeros(m,5);
end
for kc = 1:length(chains)
    chain = chains{kc};
    if isstrprop(chain(1),'upper') % chain fields start with a capital
        index_vector(1) =  entity.(chain).index;
        residues = fieldnames(entity.(chain));
        for kr = 1:length(residues) % expand over all residues
            residue = residues{kr};
            if strcmp(residue(1),'R') % these are residue fields
                index_vector(2) =  entity.(chain).(residue).index;
                % all locations are selected
                location_vector = 1:length(entity.(chain).(residue).locations);
                % rotamers overrule locations
                if length(entity.(chain).(residue).populations) > 1
                    % all rotamers are selected
                    location_vector = 1:length(entity.(chain).(residue).populations); 
                end
                atoms = fieldnames(entity.(chain).(residue));
                for ka = 1:length(atoms) % expand over all atoms
                    atom = atoms{ka};
                    if isstrprop(atom(1),'upper') % these are atom fields
                        index_vector(3) =  entity.(chain).(residue).(atom).index;
                        if entity.(chain).(residue).(atom).selected % report back
                            for kconf = 1:length(conformers)
                                index_vector(4) = conformers(kconf);
                                for kl = location_vector
                                    index_vector(5) = kl;
                                    outputs = outputs + 1;
                                    switch attribute
                                        case 'bfactor'
                                            argout{outputs} = entity.(chain).(residue).(atom).bfactor;
                                        case 'charge'
                                            argout{outputs} = entity.(chain).(residue).(atom).charge;
                                        case 'element'
                                            [~,atom_indices] = ismember(entity.index_array,index_vector,'rows');
                                            atom_index = all_atom_indices(atom_indices~=0);
                                            if ~isempty(atom_index)
                                                argout{outputs} = entity.elements(atom_index);
                                            else
                                                outputs = outputs - 1;
                                            end
                                        case 'info'
                                            info.name = atom;
                                            [~,atom_indices] = ismember(entity.index_array,index_vector,'rows');
                                            atom_index = all_atom_indices(atom_indices~=0);
                                            if ~isempty(atom_index)
                                                info.indices = index_vector;
                                                info.atom_index = atom_index;
                                                argout{outputs} = info;
                                            else
                                                outputs = outputs - 1;
                                            end
                                        case 'population'
                                            [~,atom_indices] = ismember(entity.index_array,index_vector,'rows');
                                            atom_index = all_atom_indices(atom_indices~=0);
                                            if ~isempty(atom_index)
                                                argout{outputs} = entity.populations(kconf)*double(entity.occupancies(atom_index))/100;
                                            else
                                                outputs = outputs - 1;
                                            end
                                        case 'xyz'
                                            [~,atom_indices] = ismember(entity.index_array,index_vector,'rows');
                                            atom_index = all_atom_indices(atom_indices~=0);
                                            if ~isempty(atom_index)
                                                argout{outputs} = entity.xyz(atom_index,:);
                                            else
                                                outputs = outputs - 1;
                                            end
                                        case 'ecoor'
                                            all_indices(outputs,:) = index_vector;
                                        otherwise
                                            argout = {};
                                            exceptions = {MException('get_atom:unsupported_attribute', 'Attribute %s not supported',attribute)};
                                            return
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end

if strcmp(attribute,'ecoor')
    all_indices = all_indices(1:outputs,:);
    all_indices = all_indices(1:outputs,:);
    [~,atom_indices] = ismember(entity.index_array,all_indices,'rows');
    atom_indices = all_atom_indices(atom_indices~=0);
    conformer_indices = entity.index_array(atom_indices,4);
    argout{2} = atom_indices;
    ecoor = zeros(length(atom_indices),5);
    ecoor(:,1) = entity.elements(atom_indices);
    ecoor(:,2:4) = entity.xyz(atom_indices,:);
    ecoor(:,5) = entity.populations(conformer_indices).*double(entity.occupancies(atom_indices).')/100;    
    argout{1} = ecoor;
else
    argout = argout(1:outputs);
end



