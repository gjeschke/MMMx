function [argout,exceptions] = get_residue(entity,attribute,address)
%
% GET_RESIDUE Retrieves attributes of a residue
%
%   [argout,exceptions] = GET_RESIDUE(entity,attribute)
%   Provides attribute values and possibly exceptions for residues selected
%   in an entity
%
%   [argout,exceptions] = GET_RESIDUE(entity,attribute,address)
%   Provides attribute values and possibly exceptions for residues selected
%   by an address
%
% INPUT
% entity       entity in an MMMx format, must be provided
% address      MMMx address, 'selected' refers to the current selection
%              defaults to 'selected'
% attribute    residue attribute to be retrieved, string, defaults to
%              'info'
%              attribute   output                           Type
%              ------------------------------------------------------------
%              dssp        DSSP secondary structure         char
%              info        .number      residue number      int
%                          .tlc         three-letter code   string
%              pointcoor   CA (aa) or C4' (nt) coordinate   (1,3) double
%              populations rotamer populations              (R,1) double
%              sheet       DSSP sheet information           (1,2) double
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
                if entity.(chain).(residue).selected % report back
                    outputs = outputs + 1;
                    switch attribute
                        case 'dssp'
                            if isfield(entity.(chain).(residue),'dssp')
                                argout{outputs} = entity.(chain).(residue).dssp;
                            else
                                argout{outputs} = '';
                            end
                        case 'info'
                            info.number = str2double(residue(2:end));
                            info.tlc = entity.(chain).(residue).name;
                            argout{outputs} = info;
                        case 'populations'
                            argout{outputs} = entity.(chain).(residue).populations;
                        case 'pointcoor'
                            at_index = [];
                            % check whether there is a CA atom, if so, assign
                            % its index, if not, check for C4' atom
                            if isfield(entity.(chain).(residue),'CA')
                                at_index = entity.(chain).(residue).CA.index;
                            elseif isfield(entity.(chain).(residue),'C4_') % this is how C4' is coded
                                at_index = entity.(chain).(residue).C4_.index;
                            end
                            outputs = outputs - 1;
                            if ~isempty(at_index)
                                index_vector(3) = at_index;
                                index_vector(5) = 1;
                                for kconf = 1: length(conformers)
                                    outputs = outputs + 1;
                                    index_vector(4) = conformers(kconf);
                                    [~,atom_indices] = ismember(entity.index_array,index_vector,'rows');
                                    atom_index = all_atom_indices(atom_indices~=0);
                                    argout{outputs} = entity.xyz(atom_index,:);
                                end
                            end
                        case 'pointindices'
                            at_index = [];
                            % check whether there is a CA atom, if so, assign
                            % its index, if not, check for C4' atom
                            if isfield(entity.(chain).(residue),'CA')
                                at_index = entity.(chain).(residue).CA.index;
                            elseif isfield(entity.(chain).(residue),'C4_') % this is how C4' is coded
                                at_index = entity.(chain).(residue).C4_.index;
                            end
                            outputs = outputs - 1;
                            if ~isempty(at_index)
                                index_vector(3) = at_index;
                                index_vector(5) = 1;
                                for kconf = 1: length(conformers)
                                    outputs = outputs + 1;
                                    index_vector(4) = conformers(kconf);
                                    [~,atom_indices] = ismember(entity.index_array,index_vector,'rows');
                                    atom_index = all_atom_indices(atom_indices~=0);
                                    argout{outputs} = atom_index;
                                end
                            end
                        case 'sheet'
                            if isfield(entity.(chain).(residue),'sheet')
                                argout{outputs} = entity.(chain).(residue).sheet;
                            else
                                argout{outputs} = [];
                            end
                        otherwise
                            argout = {};
                            exceptions = {MException('get_residue:unsupported_attribute', 'Attribute %s not supported',attribute)};
                            return
                    end
                end
            end
        end
    end
end

argout = argout(1:outputs);

profile viewer
