function entity = add_residue(entity,chain,resnum,residue,xyz)
%
% ADD_RESIDUE Add a residue to an entity
%
%   entity = ADD_RESIDUE(entity,residue,xyz)
%
% INPUT
% entity    MMMx entity to which the residue is to be added
% chain     chain identifier
% resnum    residue number
% residue   residue field from another entity
% xyz       coordinate array of the other entity
%
% OUTPUT
% entity    the entity with the new residue added
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

maxat = 1000; % maximum number of atoms in the residue

% extract old coordinate array of existing entity
xyz0 = entity.xyz;
occ0 = entity.occupancies;
[atnum,~] = size(xyz0); % number of atoms in the old entity

% add residue with original coordinate table entries
resname = sprintf('R%i',resnum);
entity.(chain).(resname) = residue;

xyz0 = [xyz0; zeros(maxat,3)]; % enlarge coordinate array
occ0 = [occ0; 100*ones(maxat,1,'uint8')];
atoms = fieldnames(residue);
for ka = 1:length(atoms) % loop over all atoms
    atom = atoms{ka};
    if isstrprop(atom(1),'upper') % these are atom fields
        tab_indices = residue.(atom).tab_indices;
        nat = length(tab_indices);
        coor = xyz(tab_indices,:);
        xyz0(atnum+1:atnum+nat,:) = coor;
        entity.(chain).(resname).(atom).tab_indices = atnum+1:atnum+nat;
        atnum = atnum + nat;
    end
end

entity.xyz = xyz0(1:atnum,:);
entity.occupancies = occ0(1:atnum,:);