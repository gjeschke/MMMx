function [pair_rmsd,ordering,cluster_assignment,cluster_sizes,cluster_pop,fuzziness] = cluster_sorting(pair_rmsd,pop,depth,ordering)
% [pair_rmsd,ordering,cluster_assignment,cluster_sizes,cluster_pop] = 
%              CLUSTER_SORTING(pair_rmsd,pop,ordering,options)
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
% depth         current depth, defaults to 1
% ordering      current ordering, defaults to ascending order
%
% Output:
%
% pair_rmsd             newly sorted pair_rmsd matrix
% ordering              new ordering of the structure/conformers
% cluster_assignment    assignment of conformers in original matrix to
%                       clusters
% cluster_sizes         sizes of the clusters (number of conformers)
% cluster_pop           population of the clusters
% fuzziness             if this output is requested, quality of the
%                       clustering is assessed by a fuzziness parameter,
%                       normalized intra-cluster MSD over normalized
%                       inter-cluster MSD, the lower this is, the better
%                       the clustering
% 

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% define default output
C = length(pop);
[C1,C2] = size(pair_rmsd);
if ~exist('ordering','var') || isempty(ordering)
    ordering = 1:C;
end
cluster_assignment = [];
cluster_sizes = [];
cluster_pop = [];

if C ~= C1 || C ~= C2
    return
end

if ~exist('depth','var') || isempty(depth)
    depth = 1;
end

% hierarchical clustering
pairdist = squareform(pair_rmsd);
z = linkage(pairdist);
% decrease inconsistency coefficient threshold until more than one cluster
% is found
inconsistency = 1.51;
c = ones(C,1);
while max(c) < 2 && inconsistency > 0.01
    inconsistency = inconsistency - 0.01;
    c = cluster(z,'cutoff',inconsistency);
end

% Compute populations of the clusters
pop_ordered = pop(ordering);
cluster_pop = zeros(1,max(c));
cluster_sizes = zeros(1,max(c));
for k = 1:length(c)
    cluster_pop(c(k)) = cluster_pop(c(k)) + pop_ordered(k);
    cluster_sizes(c(k)) = cluster_sizes(c(k)) + 1;
end
% Determine the cluster with maximum population
[~,poi] = max(cluster_pop);

% Find standard deviation of all clusters from the one with maximum
% population
cluster_std_dev = zeros(1,max(c));
n1 = cluster_sizes(poi);
for k = 1:max(c)
    if k ~= poi
        n2 = cluster_sizes(k);
        cluster_pair_rmsd = zeros(n1,n2);
        pop1 = zeros(1,n1);
        pop2 = zeros(1,n2);
        k1 = 0;
        k2 = 0;
        for kc1 = 1:length(c)-1
            for kc2 = kc1+1:length(c)
                if c(kc1) == poi && c(kc2) == k
                    k1 = k1 + 1;
                    pop1(k1) = pop_ordered(kc1);
                    k2 = k2 + 1;
                    pop2(k2) = pop_ordered(kc2);
                    cluster_pair_rmsd(k1,k2) = pair_rmsd(kc1,kc2);
                end
                if c(kc2) == poi && c(kc1) == k
                    k1 = k1 + 1;
                    pop1(k1) = pop_ordered(kc2);
                    k2 = k2 + 1;
                    pop2(k2) = pop_ordered(kc1);
                    cluster_pair_rmsd(k1,k2) = pair_rmsd(kc1,kc2);
                end
            end
        end
        cluster_std_dev(k) = ensemble_std(cluster_pair_rmsd,pop1,pop2);
    end
end

% Sort clusters by standard deviation form the largest (by population) cluster
[~,cluster_ordering] = sort(cluster_std_dev);
% Make the new index vector for conformers
new_ordering = zeros(1,length(ordering));
spoi = 0;
for k = 1:length(cluster_ordering) % all clusters in the correct sequence
    for kc = 1:length(c) % all conformers
        if c(kc) == cluster_ordering(k) % conformer is assigned to current cluster
            spoi = spoi + 1;
            new_ordering(spoi) = kc;
        end
    end
end

% rearrange pair rmsd matrix and population vector
pair_rmsd = pair_rmsd(new_ordering,new_ordering);
ordering = ordering(new_ordering);
cluster_assignment = c(new_ordering);
cluster_sizes = cluster_sizes(cluster_ordering);
cluster_pop = cluster_pop(cluster_ordering);
block_pointer = 0;
new_ordering = 1:length(ordering);
depth = depth + 1;
for k = 1:length(cluster_sizes)
    if cluster_sizes(k) > 2 && max(abs(c-mean(c))) > 0
        block_pair_rmsd = pair_rmsd(block_pointer+1:block_pointer+cluster_sizes(k),...
            block_pointer+1:block_pointer+cluster_sizes(k));
        block_pop = pop(ordering(block_pointer+1:block_pointer+cluster_sizes(k)));
        block_ordering = 1:cluster_sizes(k);
        [~,new_block_ordering] = cluster_sorting(block_pair_rmsd,block_pop,depth,block_ordering);
        new_ordering(block_pointer+1:block_pointer+cluster_sizes(k)) = block_pointer + new_block_ordering;
    end
    block_pointer = block_pointer+ cluster_sizes(k);
end

if nargout > 5
    n_intra = 0;
    n_inter = 0;
    intra = 0;
    inter = 0;
    for c1 = 1:C-1
        for c2 = c1+1:C
            if cluster_assignment(c1) == cluster_assignment(c2)
                n_intra = n_intra + 1;
                intra = intra + pair_rmsd(c1,c2)^2;
            else
                n_inter = n_inter + 1;
                inter = inter + pair_rmsd(c1,c2)^2;
            end
        end
    end
    intra = intra/n_intra;
    inter = inter/n_inter;
    fuzziness = intra/inter; 
end

pair_rmsd = pair_rmsd(new_ordering,new_ordering);
ordering = ordering(new_ordering);
