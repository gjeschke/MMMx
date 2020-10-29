function Rg = gyration_radius(coor)
%
% GYRATION_RADIUS Compute radius of gyration for a coordinate array
%
%   Rg = GYRATION_RADIUS(coor)
%   Radius of gyration for a molecule or polymer backbone
%
%
% INPUT
% coor  (N,3) Cartesian coordinate array
%
% OUTPUT
% Rg    radius of gyration
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

[N,~] = size(coor); % number of atoms
rc = mean(coor); % center coordinate
rel_coor = coor - repmat(rc,N,1);
Rg = sqrt(sum(sum(rel_coor.^2))/N);

