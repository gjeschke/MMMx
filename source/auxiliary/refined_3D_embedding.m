function [coor,rmsd,errcode] = refined_3D_embedding(D,Rg,weights)
% [coor,rmsd,all_rmsd] = refined_3D_embedding(D,Rg)
%
% Embeds a distance matrix D in 3D space by non-classical multi-dimensional
% scaling
%
% (c) G.Jeschke, 2024

[~,C] = size(D);

if ~exist('weights','var') || isempty(weights)
    weights = ones(C,1)/C;
end

errcode = 0;

try
    [coor,stress] = mdscale(D,3,'Weights',kron(weights,weights'));
    rmsd = sqrt(stress);
catch
    errcode = 1;
    try
        [coor,rmsd] = iterative_embedding(D);
    catch
        errcode = 2;
        coor = cmdscale(D,3);
        D_check = squareform(pdist(coor));
        rmsd = sqrt(2*sum(sum(triu(D_check-D).^2))/(C*(C-1)));
    end
end

centred = true;
inertia = get_inertia_tensor(coor,centred,weights);
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

function [coor,rmsd] = iterative_embedding(D)

[~,C] = size(D);

coor = cmdscale(D,3);
D_check = squareform(pdist(coor));
violation = sum(sum(triu(D_check-D).^2));
rmsd = 2*sqrt(2*sum(sum(triu(D_check-D).^2))/(C*(C-1)));

iter=round(violation);
if iter<1000
    iter=1000;
end
if iter>50000
    iter=50000;
end

lambda = linspace(1/iter,1,iter);

all_rmsd = zeros(1,length(lambda));

for k=1:iter
    ndmat = D + (D==0);
    D_check = squareform(pdist(coor));
    ndev = (D-D_check)./ndmat;
    if sum(sum(isnan(ndev)))
        break
    end
    coor1 = coor;
    for p = 1:C
        cvec=zeros(1,3);
        for pp=1:C
            cvec = cvec + ndev(p,pp)*(coor(p,:)-coor(pp,:));
        end
        coor1(p,:) = coor(p,:) + lambda(k)*cvec;
    end
    D_check = squareform(pdist(coor1));
    all_rmsd(k) = 2*sqrt(2*sum(sum(triu(D_check-D).^2))/(C*(C-1)));
    if all_rmsd(k) < rmsd
        rmsd = all_rmsd(k);
        coor = coor1;
    else
        break
    end
end
