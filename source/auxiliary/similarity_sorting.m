function [pair_rmsd,ordering] = similarity_sorting(pair_rmsd,pop)
% [pair_rmsd,ordering] = 
%              SIMILARITY_SORTING(pair_rmsd,pop)
%
% Sorts an ensemble of structure models/conformers according to similarity
% The central structure (with lowest mean square deviation from all others)
% tends to occur early in the sorting
% based on (iterative) hierachical clustering and ordering of clusters by
% their population-weighted mean square deviation from the larges cluster
% by population
%
% pair_rmsd     matrix of pairwise rmsd between structures (or of some
%               other distance metric), must be square and at least (3,3)
% pop           populations
%
% Output:
%
% pair_rmsd             newly sorted pair_rmsd matrix
% ordering              new ordering of the structure/conformers
% 

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

% define default output
C = length(pop);
[C1,C2] = size(pair_rmsd);

if C ~= C1 || C ~= C2
    return
end

ordering = zeros(1,length(pop));
[~,indices] = sort(pop,'descend');
current = indices(1);
ordering(1) = current;
for k = 2:length(pop)
    [~,indices] = sort(pair_rmsd(current,:));
    k2 = 1;
    while min(abs(ordering(1:k-1)-indices(k2))) == 0
        k2 = k2 + 1;
    end
    ordering(k) = indices(k2);
    current = indices(k2);
end

% rearrange pair rmsd matrix and population vector
pair_rmsd = pair_rmsd(ordering,ordering);
