function entity = asphericity(entity,selection)
%
% ASPHERICITY Computes asphericity and radius of gyration for an ensemble
%             and stores it in the entity
%
%   entity = asphericity(entity)
%   asphericity and Rg with respect to all atoms of the entity 
%
%   entity = asphericity(entity,selection)
%   asphericity and Rg with respect to a selection of atoms from the entity
%
% INPUT
% entity    entity in an MMMx format
% selection MMMx address of a selected chain or residue range
%
% OUTPUT
% entity    entity with additional fields for asphericity and Rg, the
%           fields are
%           .Rg_c           vector of radii of gyration of all conformers
%           .asphericity_c  vector of asphericity of all conformers
%           .asphericity    ensemble asphericity
%
% asphericity was defined in J Rudnick and G Gaspari 1986 J. Phys. A: Math.
% Gen. 19 L191, DOI: 10.1088/0305-4470/19/4/004
% here, we use the protein ensemble definition in
% J. Phys. Chem. B 2015, 119, 41, 12990–13001, DOI: 10.1021/acs.jpcb.5b07124

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

chemistry = load('element_attributes.mat');

Rg_all = zeros(length(entity.populations),1);
asph_all = Rg_all;
I1_squared_mean = 0;
I2_mean = 0;
for c = 1:length(entity.populations)
    % select the coordinates for which the inertia tensor is computed and the
    % corresponding elements
    if exist('selection','var') && ~isempty(selection)
        address = sprintf('{%i}%s',c,selection);
        overwrite = true;
        entity = select(entity,address,overwrite);
        [atom_indices] = get_selection(entity);
    else
        argout = get_conformer(entity,'coor',c);
        atom_indices = argout{1}.indices;
    end
    coor = entity.xyz(atom_indices,:);
    if max(atom_indices) > length(entity.elements)
        elements = entity.elements(atom_indices - (c-1)*length(entity.elements));
    else
        elements = entity.elements(atom_indices);
    end
    % compute atom masses
    masses = chemistry.pse.mass(elements).';
    % get the inertia tensor in the original frame
    centred = true;
    inertia = get_inertia_tensor(coor,centred,masses);
    I = eig(inertia);
    I1 = sum(I);
    I1_squared_mean = I1_squared_mean + entity.populations(c)*I1^2;
    I2 = I(1)*I(2) + I(1)*I(3) + I(2)*I(3);
    I2_mean = I2_mean + entity.populations(c)*I2;
    asph_all(c) = 1 - 3*I2/I1^2;   
    Rg_all(c) = gyration_radius(coor);
end
entity.rg_c = Rg_all;
entity.asphericity_c = asph_all;
entity.asphericity = 1 - 3*I2_mean/I1_squared_mean;





