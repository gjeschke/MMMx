function [pair_rmsd,ordering,cluster_assignment,cluster_sizes,cluster_pop] = cluster_analysis(pair_rmsd,pop,ordering,options)
% [pair_rmsd,ordering,cluster_assignment,cluster_sizes,cluster_pop] = 
%              cluster_analysis(pair_rmsd,pop,ordering,options)
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
% pop           populations, can be order differently from pair_rmsd
% ordering      ordering of the original structures in the pair_rmsd matrix
%               defaults to ascending order
% options       struct of options
%               .iterative  (default: true) flag indicating iterative mode
%               .verbose    (default:false) flag for verbose output to
%                           command window
%               .fname      optional file name for assignment file
%               .depth      recursion depth, defaults to 1
%               .previous   cluster number in previous iteration, defaults
%                           to 0
%
% Output:
%
% pair_rmsd             newly sorted pair_rmsd matrix
% ordering              new ordering of the structure/conformers
% cluster_assignment    assignment of conformers in original matrix to
%                       clusters
% cluster_sizes         sizes of the clusters (number of conformers)
% cluster_pop           population of the clusters
% 
% G. Jeschke, 22.12.2019

% initalize empty output
cluster_assignment = [];
cluster_sizes = [];
cluster_pop = [];

[C,n2] = size(pair_rmsd);
if C ~= n2 || C < 3
    return
end
if ~exist('ordering','var') || isempty(ordering)
    ordering = 1:C;
end

if ~exist('options','var') || isempty(options) || ~isfield(options,'iterative') || isempty(options.iterative)
    options.iterative = true;
end
if ~isfield(options,'depth') || isempty(options.depth)
    options.depth = 1;
end
if ~isfield(options,'previous') || isempty(options.previous)
    options.previous = 0;
end
if ~isfield(options,'verbose') || isempty(options.verbose)
    options.verbose = false;
end

% hierarchical clustering
pairdist = squareform(pair_rmsd);
z = linkage(pairdist);
% decrease inconsistency coefficient threshold until more than one cluster
% is found
inconsistency = 1.51;
c = ones(length(pop),1);
while max(c) < 2 && inconsistency > 0.01
    inconsistency = inconsistency - 0.01;
    c = cluster(z,'cutoff',inconsistency);
end
if max(c) < 2
    options.iterative = false;
end
if options.verbose
    fprintf(1,'%i clusters detected at inconsistency coefficient threshold %4.2f\n',max(c),inconsistency);
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
[ma,poi] = max(cluster_pop);
if options.verbose
    fprintf(1,'Cluster %i with size %i has maximum population %5.3f\n',poi,cluster_sizes(poi),ma);
end
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
        cluster_std_dev(k) = ensemble_standard_deviation(cluster_pair_rmsd,pop1,pop2);
        if options.verbose
            fprintf(1,'Cluster %i with size %i has standard deviation of %5.3f A to largest cluster %i\n',k,n2,cluster_std_dev(k),poi);
        end
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
if options.iterative
    new_ordering = 1:length(ordering);
    options.depth = options.depth + 1;
    for k = 1:length(cluster_sizes)
        options.previous = k;
        if cluster_sizes(k) > 2
            block_pair_rmsd = pair_rmsd(block_pointer+1:block_pointer+cluster_sizes(k),...
                block_pointer+1:block_pointer+cluster_sizes(k));
            block_pop = pop(ordering(block_pointer+1:block_pointer+cluster_sizes(k)));
            block_ordering = 1:cluster_sizes(k);
            [~,new_block_ordering] = cluster_analysis(block_pair_rmsd,block_pop,block_ordering,options);
            new_ordering(block_pointer+1:block_pointer+cluster_sizes(k)) = block_pointer + new_block_ordering;
        end
        block_pointer = block_pointer+ cluster_sizes(k);
    end
    pair_rmsd = pair_rmsd(new_ordering,new_ordering);
    ordering = ordering(new_ordering);
end


function ens_std = ensemble_standard_deviation(pair_rmsd,pop1,pop2,ind1,ind2)
% ens_std = ensemble_standard_deviation(pair_rmsd,pop,ind1,ind2)
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