function [measures,correlations] = analyze_ensemble(backbones,pop,options)
%
% ANALYZE_ENSEMBLE Ensemble analysis with a variety of options
%
%   [measures,correlations] = ANALYZE_ENSEMBLE(backbones,pop,options)
%   Given backbone coordinates and populations, different enesemble metrics
%   and correlation matrices are computed
%
% INPUT
% backbones     backbone structure (use get_backbones_ensemble to obtain)
%               struct with fields (chain) for all chains, each chain has
%               fields
%               .type 1 for peptide, 2 for nucleic acid, for mixed chains,
%                     the type with more backbone atoms prevails
%               .mono (1,nm) double, numbers of the residues
%               .bb   (1,C) cell of (nm,3) coordinate array for backbone
%                     atoms
% pop           (1,C) population vector for C conformers, defaults to
%               uniform populations
% options       options for output to be computed, struct with fields
%               .chain_mode     Boolean, if true, analysis for individual
%                               chains, defaults to false
%               .Rg             Boolean, if true, radius of gyration is
%                               analyzed, defaults to true
%
%
%
% OUTPUT
% measures      ensemble measures, empty if no input information
%               struct with fields for all chains (chain mode) or a single
%               field .all (chain mode false), each field has subfields
%               .Rg         radius of gyration
%               .std_Rg     standard deviation of radius of gyration  
% correlations  correlation matrices, empty if no input information
%               struct with fields for all chains (chain mode) or a single
%               field .all (chain mode false), each field has subfields
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

measures = [];
correlations = [];

chains = fieldnames(backbones);
if isempty(chains)
    return
end
C = length(backbones.(chains{1}).bb);

if ~exist('pop','var') || isempty(pop)
    pop = ones(1,C)/C;
end

if ~exist('options','var') || isempty(options) || ...
        ~isfield(options,'chain_mode') || isempty(options.chain_mode)
    options.chain_mode = false;
end

if ~isfield(options,'Rg') || isempty(options.Rg)
    options.Rg = true;
end

if ~options.chain_mode
    % determine total size of coordinate array
    N = 0;
    for kc = 1:length(chains)
        bb0 = backbones.(chains{kc}).bb{1};
        [N0,~] = size(bb0);
        N = N + N0;
    end
    for c = 1:C
        % combine coordinate arrays
        bb = zeros(N,3);
        poi = 0;
        for k = 1:length(chains)
            bb0 = backbones.(chains{k}).bb{c};
            [N0,~] = size(bb0);
            bb(poi+1:poi+N0,:) = bb0;
            poi = poi + N0;
        end
        % assign complete backbone
        backbones.all.bb{c} = bb;
    end
    % remove chain backbones
    for kc = 1:length(chains)
        backbones = rmfield(backbones,chains{kc});
    end
    chains = fieldnames(backbones);
end

