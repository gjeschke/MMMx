function entity = inertia_frame(entity,selection)
%
% INERTIA_FRAME Transforms entity into an inertia frame
%
%   entity = INERTIA_FRAME(entity)
%   inertia frame with respect to all atoms of the entity 
%
%   entity = INERTIA_FRAME(entity,selection)
%   inertia frame with respect to a selection of atoms from the entity
%
% INPUT
% entity    entity in an MMMx format
% selection MMMx address of a selected chain or residue range
%
% OUTPUT
% entity    entity with transformed coordinates
%
% After the transformation, the frame origin is the centre of gravity of 
% the selected atoms (or all atoms)
% the z axis corresponds to the maximum principal value of the inertia
% tensor, the x axis to the minimum principal value
% for a single chain, the N terminus or 5'end (first atom of the chain) 
% has a lower x and z coordinate than the C terminus or 3' end
% for multi-chain entities coordinate ordering is consistent between
% conformers (first atom has lower x and z coordinates than last atom), but
% does not have a simple physical interpretation
% for multiple conformers, each conformer is put into its own inertia frame

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

chemistry = load('element_attributes.mat');

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
        if isempty(atom_indices) && length(entity.populations) == 1
            [na,~] = size(entity.xyz);
            atom_indices = 1:na;
        end
    end
    coor = entity.xyz(atom_indices,:);
    elements = entity.elements(atom_indices);
    % compute atom masses
    masses = chemistry.pse.mass(elements).';
    % get the inertia tensor in the original frame
    centred = true;
    inertia = get_inertia_tensor(coor,centred,masses);
    % diagonalize the inertia tensor
    [evec,D] = eig(inertia);
    % sort eigenvalues in ascending order
    [~,indices] = sort(diag(D));
    % sort eigenvectors accordingly
    evecs = evec(:,indices);
    % determine the centre of gravity of the selected atoms
    xyzc = mean(masses.*coor)/mean(masses);
    
    % determine indices of all atoms of this conformer, if there was a
    % selection
    if ~exist('selection','var') || isempty(selection)
        argout = get_conformer(entity,'coor',c);
        atom_indices = argout{1}.indices;
        if isempty(atom_indices) && length(entity.populations) == 1
            [na,~] = size(entity.xyz);
            atom_indices = 1:na;
        end
    end
    % make the transformation
    xyz = entity.xyz(atom_indices,:);
    [n,~] = size(xyz); % number of atoms
    % put orgin at centre of gravity of the selected atoms
    xyz = xyz - repmat(xyzc,n,1);
    % transform into inertia frame
    xyz = xyz*evecs;
    % invert coordinates so that N terminus has lower z and x coordinates
    if xyz(1,3) > xyz(end,3) % the z coordinate must be inverted
        xyz(:,3) = -xyz(:,3);
        % this requires that either the x or the y coordinate is also inverted
        % (keep chirality)
        if xyz(1,1) > xyz(end,1) % in this case, the x coordinate should be inverted
            xyz(:,1) = -xyz(:,1);
        else % otherwise, the y coordinate is inverted
            xyz(:,2) = -xyz(:,2);
        end
    else % the z coordinate was not inverted
        % if the x coordinate must be inverted, the y coordinate must be
        % inverted, too, to keep chirality
        if xyz(1,1) > xyz(end,1) % in this case, the x coordinate should be inverted
            xyz(:,1) = -xyz(:,1);
            xyz(:,2) = -xyz(:,2);
        end
    end
    entity.xyz(atom_indices,:) = xyz;
end





