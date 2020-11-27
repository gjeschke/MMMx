function Rh = hydrodynamic_radius(coor,N)
%
% HYDRODYNAMIC_RADIUS Compute hydrodynamic radius for a coordinate array
%
%   Rh = HYDRODYNAMIC_RADIUS(coor)
%   Hydrodynamic radius for a molecule or polymer backbone
%
%   based on: M. Nygaard, B.B. Kragelund, E. Papaleo, K. Lindorff-Larsen,
%   Biophys. J. (2017) 113, 550-557.
%
% INPUT
% coor  (M,3) Cartesian coordinate array
% N     number of residues
%
% OUTPUT
% Rh    hydrodynamic radius
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% Nygaard parameters

alpha1 = 0.216;
alpha2 = 4.06;
alpha3 = 0.821;

% radius of gyration
[M,~] = size(coor); % number of atoms
rc = mean(coor); % center coordinate
rel_coor = coor - repmat(rc,M,1);
Rg = sqrt(sum(sum(rel_coor.^2))/M);

% Eq. (7) of the Nygaard paper solved for Rh
Rh = (alpha3-Rg)*(sqrt(N)-N^(1/3))/(alpha1*(alpha2*N^(1/3)-Rg));

