function energy = energy_LJ(coor1,coor2,options)
%
% ENERGY_LJ computes Lennard-Jones energy between two sets of extended
%           coordinates
%
%   energy = ENERGY_LJ(coor1,coor2,options)
%   Given two sets of coordinates with associated atomic numbers and at
%   least Lennard-Jones definitions of a force field in options, computes
%   Lennard-Jones interaction energy in J/mol
%
% INPUT
% coor1    (N,4) double for N atoms, first column atomic numbers, columns
%           2-4 Cartesian coordinates
% coor2    (M,4) double for M atoms, first column atomic numbers, columns
%           2-4 Cartesian coordinates
% options  struct with run options, force field is obligatory
%          .ff       force field
%          .ff.LJ_r  Lennard-Jones (van-der-Waals) radii by atomic number
%          .ff.D     Lennard-Jones energies by atomic number
%          .f_factor scaling factors, (1,2) double
%                    .f_factor(1) scales van-der-Waals radius, typically 0.5 
%                    .f_factor(2) scales attraction w.r.t. repulsion term,
%                                 typically 1, 2 for "sticky" labels
%
% OUTPUT
% energy    interaction energy in J/mol

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2010-2020: Yevhen O. Polyhach, Stefan Stoll, Gunnar Jeschke


conv_factor=(4.185*1000); % conversion from kcal/mol to J/mol

% set default f-factor and no enhancement of attractive term
if ~isfield(options,'f_factor')
    options.f_factor = [0.5,1];
end
if ~isfield(options,'cutoff')
    options.cutoff = 10;
end

[N,~] = size(coor1); % get sizes of the coordinates arrays
[M,~] = size(coor2);

a = coor1(:,2:4);
b = coor2(:,2:4);
a2 = repmat(sum(a.^2,2),1,M);
b2 = repmat(sum(b.^2,2),1,N).';
pair_dist = sqrt(abs(a2 + b2 - 2*a*b.'));
pair_dist = reshape(pair_dist,1,M*N);

Rmin2_i0 = options.ff.LJ_r(coor1(:,1));
eps_i = options.ff.LJ_D(coor1(:,1))*conv_factor;
Rmin2_j0 = options.ff.LJ_r(coor2(:,1));
eps_j = options.ff.LJ_D(coor2(:,1))*conv_factor;

eps = sqrt(eps_i(:)*eps_j);
eps = reshape(eps,1,M*N);
Rmin2_i = Rmin2_i0*options.f_factor(1); % soften LJ with the forgive factor
Rmin2_j = Rmin2_j0*options.f_factor(1); % soften LJ with the forgive factor
Rmin = sqrt(Rmin2_i(:)*Rmin2_j); % UFF has geometric average combination rule
eps = eps(pair_dist<=options.cutoff);
Rmin = Rmin(pair_dist<=options.cutoff);
pair_dist = pair_dist(pair_dist<=options.cutoff);
q = (Rmin./pair_dist).^2;
q = q.*q.*q;
energies = eps.*(q.^2-2*q*options.f_factor(2));
energy = sum(energies);

