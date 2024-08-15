function [ensemble,info,assignment] = cluster_ensemble(entity,options)
%
% CLUSTER_ENSEMBLE Reduce ensemble size by clustering
%
%   entity = CLUSTER_ENSEMBLE(entity)
%   Reduces ensemble size by clustering based on distance root mean square
%   deviation between conformers
%
% INPUT
% entity        input entity, should have between 10 and 50000 conformers
% options       options for ensemble reduction
%               .factor         reduction factor, defaults to 5, output
%                               ensemble is five times smaller
%               .size           size of the output ensemble, overwrites
%                               .factor if present and not empty, no
%                               default
%               .chain          restricts analysis to a single chain,
%                               defaults to all chains
%               .range          restricts analysis to a residue range,
%                               range is interpreted as a list of residues 
%                               if its elements are negative 
% OUTPUT
% ensemble      [n,2] array, first column a conformer indices to the
%               original ensemble, second column new populations
% info          information on old and new ensemble
%               .entropies  [2,1] Shannon entropies of old and new ensemble
%               .widths     [2,1] widths of old and new ensemble
%               .similarity similarity measure between ensembles
%               .resolution maximum distance rmsd between conformers that
%                           were assigned to the same cluster
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

if ~exist('options','var') || isempty(options) || ~isfield(options,'factor')
    options.factor = 5;
end

if ~isfield(options,'size') || isempty(options.size)
    options.size = ceil(length(entity.populations)/options.factor);
end
pop = entity.populations;
info.entropies = [-sum(pop.*log10(pop/sum(pop)))/sum(pop),0];

dmat = pair_drms_matrix(entity);
var11 = kron(entity.populations,entity.populations').*dmat.^2;
width1 = sqrt(sum(sum(var11)));
info.widths = [width1,0];

% compute linkage
Z = linkage(dmat,'average');
% cluster
assignment = cluster(Z,'maxclust',options.size);
% select representative model for each cluster
selected_models = zeros(options.size,1);
all_indices = 1:length(pop);
cluster_resolution = 0;
pop2 = zeros(options.size,1);
for km = 1:options.size
    % determine indices corresponding to cluster km
    indices = all_indices(assignment==km);
    pop2(km) = sum(pop(indices));
    % extract distance matrix for this cluster
    dmat_cluster = dmat(indices,indices);
    % compute mean square sum for each model in cluster with respect to
    % all other models
    msq_models = sum(dmat_cluster.^2);
    % find minimum mean square deviation and index of this model in
    % cluster
    [~,index] = min(msq_models);
    selected_models(km) = indices(index);
    % maximum rooot mean square distance of another model in this
    % cluster from selected model
    max_drmsd = max(dmat_cluster(index,:));
    if max_drmsd > cluster_resolution
        cluster_resolution = max_drmsd;
    end
end


dmat_ensemble2 = dmat(selected_models,selected_models);
var22 = kron(pop2,pop2').*dmat_ensemble2.^2;
width2 = sqrt(sum(sum(var22)));
info.widths(2) = width2;
info.entropies(2) = -sum(pop2.*log10(pop2/sum(pop2)))/sum(pop2);
overlap = dmat(selected_models,:);
var12 = kron(pop2,pop').*overlap.^2;
info.similarity = sqrt((width1*width2)/sum(sum(var12)));
info.resolution = cluster_resolution;

ensemble = [selected_models pop2];
