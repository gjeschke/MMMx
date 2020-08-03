function [argout,exceptions] = get_atom(entity,attribute,address)
%
% GET_ATOM Retrieves attributes of an atom
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
%              'info', all loactions of the atom are included
%              attribute   output                           Type
%              ------------------------------------------------------------
%              bfactor     B factor                         double
%              charge      atom charge                      int
%              coor        Cartesian coordinate array       (N,3) double
%              element     atomic number                    int8
%              info        .name          atom name         string
%                                         R# for a rotamer  # is a number
%                          .indices       index vector      (1,5) uint16
%                          .atom_index    atom array index  int
%              population  rotamer population or occupancy  double
%              xyz         Cartesian coordinates            (1,3) double
%
% OUTPUT
% argout       cell array of outputs, see above
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

profile on

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
                                    [~,atom_indices] = ismember(entity.index_array,index_vector,'rows');                                    
                                    atom_index = all_atom_indices(atom_indices~=0);
                                    if ~isempty(entity.xyz(atom_index,:))
                                        outputs = outputs + 1;
                                        switch attribute
                                            case 'bfactor'
                                                argout{outputs} = entity.(chain).(residue).(atom).bfactor;
                                            case 'charge'
                                                argout{outputs} = entity.(chain).(residue).(atom).charge;
                                            case 'element'
                                                argout{outputs} = entity.elements(atom_index);
                                            case 'info'
                                                info.name = atom;
                                                info.indices = index_vector;
                                                info.atom_index = atom_index;
                                                argout{outputs} = info;
                                            case 'population'
                                                argout{outputs} = entity.populations(kconf)*double(entity.occupancies(atom_index))/100;
                                            case {'xyz','coor'}
                                                argout{outputs} = entity.xyz(atom_index,:);
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
end

argout = argout(1:outputs);

if strcmp(attribute,'coor')
    coor = cell2mat(argout);
    coor = reshape(coor,3,outputs).';
    argout = {coor};
end

profile viewer
