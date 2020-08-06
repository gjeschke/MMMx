function [entity,exceptions] = set_location(entity,attribute,argin,address)
%
% SET_LOCATION Assigns attributes to atom locations 
%
%   [entity,exceptions] = SET_LOCATION(entity,attribute,argin)
%   Modifies an attribute of selected locations in entity to the values 
%   provided in argin
%
%   [entity,exceptions] = SET_LOCATION(entity,attribute,argin,address)
%   Modifies an attribute of the locations in entity specified by address
%   to the values provided in argin
%
% INPUT
% entity       entity in an MMMx format, must be provided
% attribute    location attribute to be modified, string, must be
%              provided
%              attribute   output                           Type
%              ------------------------------------------------------------
%              element     atomic number                    int8
%              population  rotamer population or occupancy  double
%              xyz         Cartesian coordinates            (1,3) double
% argin        new values for the attribute, must be provided, cell array
%              of inputs, see table above
% address      MMMx address, 'selected' refers to the current selection
%              defaults to 'selected'
%
% OUTPUT
% entity       the input entity with  modified attributes
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
exceptions = {[]};

% set default arguments
if ~exist('address','var')
    address = 'selected';
end
if ~exist('attribute','var') || isempty(attribute)
    exceptions = {MException('set_location:missing_attribute', 'Attribute must be provided')};
end
if ~exist('argin','var')
    exceptions = {MException('set_location:missing_input', 'Input arguments must be provided')};
end

% select the objects by the provided address
entity = select(entity,address,true);

inputs = 0; % counter for the number of outputs
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
                atoms = fieldnames(entity.(chain).(residue));
                for ka = 1:length(atoms) % expand over all atoms
                    atom = atoms{ka};
                    if isstrprop(atom(1),'upper') % these are atom fields
                        index_vector(3) =  entity.(chain).(residue).(atom).index;
                        if entity.(chain).(residue).(atom).selected % report back
                            for kconf = 1:length(conformers)
                                index_vector(4) = conformers(kconf);
                                for kl = 1:length(entity.(chain).(residue).(atom).selected_locations)
                                    index_vector(5) = entity.(chain).(residue).(atom).selected_locations(kl);
                                    [~,atom_indices] = ismember(entity.index_array,index_vector,'rows');                                    
                                    atom_index = all_atom_indices(atom_indices~=0);
                                    if ~isempty(entity.xyz(atom_index,:))
                                        inputs = inputs + 1;
                                        switch attribute
                                            case 'element'
                                                entity.elements(atom_index) = argin{inputs};
                                            case 'population'
                                                entity.occupancies(atom_index) =  uint8(100*argin{inputs}/entity.populations(kconf));
                                            case {'xyz'}
                                                entity.xyz(atom_index,:) = argin{inputs};
                                            otherwise
                                                exceptions = {MException('set_location:unsupported_attribute', 'Attribute %s not supported',attribute)};
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
