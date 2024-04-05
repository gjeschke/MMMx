function [pair_drmsd,Rg,ordering] = similarity_sorting(pair_drmsd,Rg)
% [pair_drmsd,Rg,ordering] = 
%              SIMILARITY_SORTING(pair_drmsd,Rg)
%
% Sorts an ensemble of structure models/conformers according to similarity
% The central structure (with lowest mean square deviation from all others)
% tends to occur early in the sorting
% based on (iterative) hierachical clustering and ordering of clusters by
% their population-weighted mean square deviation from the larges cluster
% by population
%
% pair_drmsd    matrix of pairwise distance rmsd between structures (or of
%               some other distance metric), must be square and at least (3,3)
% Rg            vector of radii of gyration
%
% Output:
%
% pair_rmsd             newly sorted pair_rmsd matrix
% Rg                    newly ordered vector of radii of gyration
% ordering              new ordering of the structure/conformers
% 

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

% define default output
[C1,C2] = size(pair_drmsd);

if C1 ~= C2
    return
end

% Find the maximum distance rmsd between two conformers
max_drmsd = max(max(pair_drmsd));

% Find the indices of these conformers
[c1,c2] = find(pair_drmsd == max_drmsd);
if length(c1) > 1
    c1 = c1(1);
end
if length(c2) > 1
    c2 = c2(1);
end
if Rg(c1) > Rg(c2)
    xc = c2;
    c2 = c1;
    c1 = xc;
end

Dstart = pair_drmsd(c1,:) + 1e-12;
Dend = pair_drmsd(c2,:) + 1e-12;
ranking = Dstart./Dend;
[~,ordering] = sort(ranking);

ordered = false;
iteration = 0;
max_iterations = 10000;
all_X = zeros(1,max_iterations);
while ~ordered && iteration < max_iterations
    Xsum = 0;
    new_ordering = ordering;
    ordered = true;
    iteration = iteration + 1;
    for c = 2:C1-1
        rLA = pair_drmsd(ordering(1),ordering(c)); 
        rLB = pair_drmsd(ordering(1),ordering(c+1));
        rAH = pair_drmsd(ordering(end),ordering(c)); 
        rBH = pair_drmsd(ordering(end),ordering(c+1));
        X = rLB^2 + rAH^2 - rBH^2 - rLA^2;
        Xsum = Xsum + X;
        if X < 0
            ordered = false;
            new_ordering(c) = ordering(c+1);
            new_ordering(c+1) = ordering(c);
        end
    end
    all_X(iteration) = Xsum;
    ordering = new_ordering;
end
% rearrange pair rmsd matrix and population vector
pair_drmsd = pair_drmsd(ordering,ordering);
Rg = Rg(ordering);
