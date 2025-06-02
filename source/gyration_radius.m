function Rg = gyration_radius(coor,weights)
%
% GYRATION_RADIUS Compute radius of gyration for a coordinate array
%
%   Rg = GYRATION_RADIUS(coor)
%   Radius of gyration for a molecule or polymer backbone
%
%
% INPUT
% coor      (N,~) Cartesian coordinate array, this also works in Euclidean
%           hyperspace
% weights   vector (N,1) of weights, optional, defaults to uniform weights
%
% OUTPUT
% Rg    radius of gyration
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020-2025: Gunnar Jeschke

[N,~] = size(coor); % number of atoms

if ~exist('weights','var') || isempty(weights)
    weights = ones(1,N);
end
weights = weights/sum(weights);
[n,m] = size(weights);
if n > m
    weights = weights';
end

rc = weights*coor; % center coordinate
rel_coor = coor - repmat(rc,N,1);
Rg = sqrt(sum(weights.*sum(rel_coor'.^2)));

