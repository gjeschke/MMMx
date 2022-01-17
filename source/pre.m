function [all_Gamma2,all_pre] = pre(coor,label,td,taur,taui,R2_dia,larmor)
%
% PRE Computation of paramagnetic relaxation enhancement 
%
%   [all_Gamma2,all_pre] = PRE(coor,label_coor,taur,taui,R2_dia)
%   Given nuclear and spin label coordinates and populations, the protein
%   rotational correlatuon time and label internal correlation time and the
%   nuclear transverse relaxation rate for the diamagnetically labeled
%   protein, relaxation enhancement rates and intensity ratios I_para/I_dia 
%   are estimated
%
%   see: Tesei et al. https://doi.org/10.1101/2020.08.09.243030 
%        for description of the approach
%
% INPUT
% coor      double (N,3) atom coordinates of the nuclear spins [Angstroem]
% label     struct for spin label with fields 
%           .coor   double (R,3) coordinates for R rotamers [Angstroem]
%           .pop    double (1,R) rotamer populations
% td        double, total INEPT time in HSQC [s]
% taur      double, rotational correlation time of the protein [s]
% taui      double, correlation time of internal label motion [s]
% R2_dia    nuclear transverse relataion rate (s^-1), can be either one double
%           value applying to all nuclei or a vector (1,N) of values for
%           individual nuclei [s^{-1)]
% larmor    larmor frequency [MHz]
%
% OUTPUT
% all_Gamma2    double(N,1) of relaxation enhancement rates
% all_pre       double (N,1) of intensity ratios I_para/I_dia
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2010-2020: Yevhen Polyhach, Gunnar Jeschke

% constants

mu_B = 9.274009994e-24; % Bohr magneton
g = 2.0059; % nitroxide g value)
gamma_H = 2.675222005e8; % proton gyromagnetic ratio
mu0o4pi = 1e-7; % vacuum permeability over 4 pi
S = 1/2;
prefac = mu0o4pi*gamma_H*g*mu_B;
prefac = S*(S+1)*prefac^2/15; % Eq. (8), Tesei et al.
angular_frq = 2*pi*1e6*larmor;

taut = 1/(1/taur + 1/taui);

% determine size of the problem
[N,~] = size(coor);
R = length(label.pop);

% initialize output
all_Gamma2 = NaN(N,1);
all_pre = NaN(N,1);

for n = 1:N % loop over nuclei
    rm3_avg = 0;
    rm6_avg = 0;
    S2_angular = 0;
    for i = 1:R % loop over first rotamer
        ri_vec = 1e-10*(label.coor(i,:) - coor(n,:)); % Angstroem to m
        ri = norm(ri_vec);
        rm3_avg = rm3_avg + label.pop(i)*ri^(-3);
        rm6_avg = rm6_avg + label.pop(i)*ri^(-6);
        for j = 1:R % loop over second rotamer
            rj_vec = 1e-10*(label.coor(j,:) - coor(n,:)); % Angstroem to m
            rj = norm(rj_vec);
            ang_term = sum(ri_vec.*rj_vec)/(ri*rj);
            S2_angular = S2_angular + (3*ang_term^2-1)*label.pop(i)*label.pop(j)/2;
        end
    end
    S2_radial = rm3_avg^2/rm6_avg;
    S2 = S2_radial*S2_angular;
    J0 = S2*taur + (1-S2)*taut; % Eq. (9)  Tesei et al.
    J_larmor = S2*taur/(1+(taur*angular_frq)^2) +(1-S2)*taut/(1+(taut*angular_frq)^2); 
    Gamma2 = prefac*rm6_avg*(4*J0 + 3*J_larmor); % Eq. (8) Tesei et al.
    all_Gamma2(n) = Gamma2;
    all_pre(n) = R2_dia*exp(-td*Gamma2)/(R2_dia+Gamma2);
end