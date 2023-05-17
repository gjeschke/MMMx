function D = dunn_index(dmat,assignment)
% D = DUNN_INDEX(dmat,assignment)
%
% Computes the Dunn index of a clustering from a distance matrix and an
% assignment
%
% Input:
%
% dmat          matrix of pairwise "distances" between observations
%               (conformers)
% assignment    assignment if observations to clusters
%
% Output:
%
% D             Dunn index
% 

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

% define default output

min_sep = 1e20;
max_dia = 0;
for k = 1:max(assignment)
    cluster = find(assignment == k);
    others = find(assignment ~= k);
    for c1 = 1:length(cluster)-1
        for c2 = c1+1:length(cluster)
            if dmat(cluster(c1),cluster(c2)) > max_dia
                max_dia = dmat(cluster(c1),cluster(c2));
            end
        end
    end
    for c1 = 1:length(cluster)
        for c2 = 1:length(others)
            if dmat(cluster(c1),others(c2)) < min_sep
                min_sep = dmat(cluster(c1),others(c2));
            end
        end
    end
end

D = min_sep/max_dia;