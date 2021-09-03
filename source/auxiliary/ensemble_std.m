function ens_std = ensemble_std(pair_rmsd,pop1,pop2,ind1,ind2)
% ens_std = ensemble_std(pair_rmsd,pop,ind1,ind2)
% population-weighted standard deviation of ensembles
%
% pair_rmsd     matrix of coordinate mean-square deviations between
%               ensemble members, square matrix for a single ensemble,
%               rectangular (n1,n2) for two ensembles, if two ensembles
%               have the same number of elements and cross standard, 
%               deviation is requested, index ranges must be given
% pop1          vector (n1,1) of populations, defaults to uniform
%               populations (1/n1)
% pop2          vector (n2,1) of populations, defaults to uniform
%               populations (1/n2), if pop2 is not given or empty, and the
%               matrix pair_rmsd is square, a single ensemble is assumed
% ind1          index range along first dimension, defaults to all
%               n1 conformers
% ind2          index range along second dimension, defaults to all
%               n2 conformers

[n1,n2] = size(pair_rmsd);

if ~exist('pop1','var') || isempty(pop1)
    pop1 = ones(1,n1)/n1;
end

% determine computation mode
single_mode = false;
if ~exist('pop2','var') || isempty(pop2)
    pop2 = ones(1,n2)/n2;
    if n1 == n2
        single_mode = true;
    end
end
if ~exist('ind1','var') || isempty(ind1)
    ind1 = 1:n1;
end
if ~exist('ind2','var') || isempty(ind2)
    ind2 = 1:n2;
end

popsum = 0; % sum of populations
msqsum = 0; % sum of mean square deviations

for k1 = ind1
    for k2 = ind2
        if ~single_mode || k1 ~= k2
            msqsum = msqsum + pop1(k1)*pop2(k2)*pair_rmsd(k1,k2)^2;
            popsum = popsum + pop1(k1)*pop2(k2);
        end
    end
end

ens_std = sqrt(msqsum/popsum);