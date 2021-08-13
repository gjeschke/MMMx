function entity = transform_structure(entity,CA_model,new_CA_coor,drag_list,conformer)
% new_entity = TRANSFORM_STRUCTURE(entity,CA_model,new_CA_coor,drag_list)
%
% transforms an entity by deformation so that it fits a new Calpha trace
%
% Input:
% entity        entity (structure) in MMMx:atomic format
% CA_model      Calpha model for m residues, struct with fields
%               .coor       (m,3) double Cartesian coordinates of CA atoms
%               .masses     (m,1) double, residue masses in g/mol (unused)
%               .bfactors   (m,1) double, B factors (unused, can be NaN)
%               .chains     (m,1) double, chain indices
%               .addresses  (m,1) cell of strings with residue addresses
% new_CA_coor   (N,3) double array of new Cartesian CA coordinates
% drag_list     list of addresses of non-petide residues that should be
%               dragged along with their nearest peptide residue (optional)
% conformer     conformer to be transformed, defaults to 1
%
% Output:
%
% entity        deformed structure, entity in MMMx:atomic format

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% G. Jeschke, 2010-2021

if ~exist('conformer','var') || isempty(conformer)
    conformer = 1;
end

if ~exist('drag_list','var')
    drag_list = {};
end

[m,~] = size(new_CA_coor);
for res = 1:m
    transmat = get_residue_transformation(res,CA_model,new_CA_coor);
    [ctag,rtag] = split_address(CA_model.addresses{res});
    atoms = fieldnames(entity.(ctag).(rtag)); % list of potential atom names
    for ka = 1:length(atoms) % expand over all atoms
        atom = atoms{ka};
        if isstrprop(atom(1),'upper') % these are atom fields
            index = entity.(ctag).(rtag).(atom).tab_indices(conformer);
            % transform atom coordinate
            xyz = [entity.xyz(index,:) 1];
            xyz = transmat*xyz';
            xyz = xyz';
            % store atom coordinate
            entity.xyz(index,:) = xyz(:,1:3);
        end
    end
end

for k = 1:length(drag_list)
    if ~isempty(drag_list{k})
        [ctag,rtag] = split_address(drag_list{k});
        atoms = fieldnames(entity.(ctag).(rtag)); % list of potential atom names
        % determine the mean coordinate of this residue
        res_coor = zeros(1,3);
        n_atoms = 0;
        for ka = 1:length(atoms) % expand over all atoms
            atom = atoms{ka};
            if isstrprop(atom(1),'upper') % these are atom fields
                n_atoms = n_atoms + 1;
                index = entity.(ctag).(rtag).(atom).tab_indices(conformer);
                % add atom coordinate
                res_coor = res_coor + entity.xyz(index,:);
            end
        end
        res_coor = res_coor/n_atoms; % mean coordinate
        % find the closest peptide residue
        diff = sum((CA_model.coor - res_coor).^2,2);
        [~,res] = min(diff);
        % affine transformation of the closest peptide residue
        transmat = get_residue_transformation(res,CA_model,new_CA_coor);
        for ka = 1:length(atoms) % expand over all atoms
            atom = atoms{ka};
            if isstrprop(atom(1),'upper') % these are atom fields
                index = entity.(ctag).(rtag).(atom).tab_indices(conformer);
                % transform atom coordinate
                xyz = [entity.xyz(index,:) 1];
                xyz = transmat*xyz';
                xyz = xyz';
                % store atom coordinate
                entity.xyz(index,:) = xyz(:,1:3);
            end
        end
    end
end

function transmat = get_residue_transformation(res,CA_model,new_CA_coor)
% provides the affine transformation matrix for a residue
%
% Input:
%
% res   residue number as pointer into CA_model.coor and new_CA_coor
% CA_model      Calpha model for m residues, struct with fields
%               .coor       (m,3) double Cartesian coordinates of CA atoms
%               .masses     (m,1) double, residue masses in g/mol (unused)
%               .bfactors   (m,1) double, B factors (unused, can be NaN)
%               .chains     (m,1) double, chain indices
%               .addresses  (m,1) cell of strings with residue addresses
% new_CA_coor   (m,3) double array of new Cartesian CA coordinates
%
% Output:
%
% transmat      4 x 4 affine transformation matrix

local_template = zeros(5,3);
local_template_0 = zeros(5,3);
trans = new_CA_coor(res,:) - CA_model.coor(res,:);
% only translation is applied, if a local frame cannot be defined
transmat = affine('translation',trans); % shift to centroid of structure to be fitted

index1 = get_residue_number(CA_model.addresses{res});

[m,~] = size(new_CA_coor);

% make a local template to fit rotation and translation
poi = 0;
for kk =- 2:2
    if res + kk > 0 && res + kk <= m % is addressed residue a network point?
        index2 = get_residue_number(CA_model.addresses{res + kk});
        diff = index2 - index1;
        if diff == kk % is addressed residue part of a continuous segment?
            poi = poi+1;
            local_template(poi,:) = new_CA_coor(res + kk,:);
            local_template_0(poi,:) = CA_model.coor(res + kk,:);
        end
    end
end
if poi > 3 % found sufficient number of points to determine local rotation and translation
    [~,~,transmat] = rmsd_superimpose(local_template(1:poi,:),local_template_0(1:poi,:));
elseif poi == 3
    [~,~,transmat] = superimpose_3points(local_template(1:poi,:),local_template_0(1:poi,:));
end

function num = get_residue_number(address)

poi = strfind(address,')');
if isempty(poi)
    poi = 0;
end
num = str2double(address(poi+1:end));

function [ctag,rtag] = split_address(address)

p1 = strfind(address,'(');
p2 = strfind(address,')');
ctag = address(p1+1:p2-1);
rtag = strcat('R',address(p2+1:end));

