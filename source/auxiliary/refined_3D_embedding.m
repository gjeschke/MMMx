function [coor,rmsd] = refined_3D_embedding(D,Rg,populations)
% [coor,rmsd,all_rmsd] = refined_3D_embedding(D,Rg)
%
% Embeds a distance matrix D in 3D space by non-classical multi-dimensional
% scaling
%
% (c) G.Jeschke, 2024

[~,C] = size(D);

[coor,stress] = mdscale(D,3,'Weights',kron(populations,populations'));
rmsd = sqrt(stress);

centred = true;
inertia = get_inertia_tensor(coor,centred,populations);
% diagonalize the inertia tensor
[evec,ID] = eig(inertia);
% sort eigenvalues in ascending order
[~,indices] = sort(diag(ID));
% sort eigenvectors accordingly
evecs = evec(:,indices);
% determine the centre of gravity of the conformer set
coorc = mean(coor);
% put orgin at centre of gravity of the conformer set
coor = coor - repmat(coorc,C,1);
% transform into inertia frame
coor = coor*evecs;

[~,Rg_sorting] = sort(Rg);
xsorting = coor(Rg_sorting,1);
invertx = sum((1:C)'.*xsorting) < 0;
ysorting = coor(Rg_sorting,2);
inverty = sum((1:C)'.*ysorting) < 0;

% invert coordinates so that lowest Rg terminus has lower x and y coordinates
if invertx % the x coordinate must be inverted
    coor(:,1) = -coor(:,1);
    % this requires that either the y or the z coordinate is also inverted
    % (keep chirality)
    if inverty % in this case, the y coordinate should be inverted
        coor(:,2) = -coor(:,2);
    else % otherwise, the z coordinate is inverted
        coor(:,3) = -coor(:,3);
    end
else % the x coordinate was not inverted
    % if the y coordinate must be inverted, the z coordinate must be
    % inverted, too, to keep the frame right-handed
    if inverty % in this case, the x coordinate should be inverted
        coor(:,2) = -coor(:,2);
        coor(:,3) = -coor(:,3);
    end
end