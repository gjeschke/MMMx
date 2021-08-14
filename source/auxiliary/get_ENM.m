function ENM_param = get_ENM(parametrization,imANM,mass_weighting)
% ENM_param = GET_ENM
%
% returns parameters of the default elastic network models of MMMx
% the default elastic network model of MMM is implemented here
%
% Input:
% parametrization   optional override of the default parametrization
%                   default is 'ed-ENM', other options are 'cutoff10', 
%                   'cutoff13', 'Hinsen', 'Jeschke', 'ed-ENM-p'
% imANM             flag, if true, anisotropic membrane protein network by
%                   Lezon/Bahar is selected, defaults to false
% mass_weighting    flag, if true, the Hessian is residue-mass weighted,
%                   defaults to false, which strangely appears to work
%                   better
%
% Output:
%
% ENM_param         struct with definition of the elastic network model

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% G. Jeschke, 2010-2021

% set default parametrization
if ~exist('parametrization','var') || isempty(parametrization)
    ENM_param.parametrization = 'ed-ENM'; % current best tested choice
else
    ENM_param.parametrization = parametrization;
end

if ~exist('imANM','var') || isempty(imANM)
    ENM_param.imANM = false; % current best tested choice
else
    ENM_param.imANM = imANM;
end

if ~exist('mass_weighting','var') || isempty(mass_weighting)
    ENM_param.mass_weighting = false; % current best tested choice
else
    ENM_param.mass_weighting = mass_weighting;
end

ENM_param.gamma = 2.06; % generic force constant divided by k_B T in ?^{-2}
ENM_param.rc= 7.3; % cutoff distance in ?, should be 7.3
ENM_param.rc_ANM = sqrt(3)*ENM_param.rc; % cutoff distance for anisotropic network model
ENM_param.p_GNM = 2; % exponent for distance-dependent force constants GNM
ENM_param.p_GNM_gamma = 2.06; % exponent for distance-dependent force constants GNM
ENM_param.p_ANM = 6; %  exponent for distance-dependent force constants ANM
ENM_param.p_ANM_gamma = 2.06; %  exponent for distance-dependent force constants ANM
ENM_param.imANM_scaling = 16; % see T. R. Lezon, I. Bahar, Biophys. J. 2012, 102, 1331-1340

ENM_param.fit_basis = 20; % number of slow modes used in structure fitting
ENM_param.fix_local = 2; % number of neighbors to which distances are fixed during fitting,
                       % 1 stabilizes only C_alpha-C_alpha distances
                       % 2 stabilizes C_alpha-C_alpha distances as well as
                       % C_alpha-C_alpha-C_alpha angles
		               % -2 stabilizes only C_alpha-C_alpha-C_alpha angles,
		               % assuming that the ENM stabilizes C_alpha-C_alpha
		               % distances
                       % 3 stabilizes local secondary structure (not yet
                       % implemented)
ENM_param.tol_local = 0.05; % tolerance (r.m.s.d. in ?) for local distortions
ENM_param.fix_force = 10; % see Zheng/Brooks (2006), below Eq. (4)
ENM_param.sec_force = 0.02; % see Zheng/Brooks (2006), below Eq. (4), different 
                        % restoring force for second neighbor
ENM_param.tol_dist = 3.0; % tolerance of distance restraints in structure fitting
ENM_param.tol_model = 0.5; % expected uncertainty of the model
ENM_param.track = false;   % initial mode set is tracked during fitting, if true, 
                        % otherwise the lowest-energy modes of new diagonalization are used
ENM_param.cycles = 100; % maximum number of iterations in ENM-based fitting
ENM_param.diagonalize = false; % if true, always recompute normal modes by Hessian rediagonalization
ENM_param.reorientate = false; % if true, always use reorientation of normal modes rather than recomputation
ENM_param.update_frequency = 5; % number of cycles, after which normal modes are updated by Hessian diagonalization
ENM_param.tif = 0.075; % threshold for active space extension
ENM_param.mmax = 30; % iteration at which active space achieves maximum dimension, can be overwritten in constraint file
ENM_param.relax = false;
ENM_param.bond_max = [4.00,7.80,11.20,14.60,18.00,21.40];   % maximum CA-CA distances for setting a force
                     % at residue number differences 1... 6


kcal2kJ=4.1868; % kcal to kJ
nm2A=10; % nanometer to ?
RT=8.314e-3*298; % 8.314 J mol^(-1) K^(-1)? 298 K
rCaCa=3.8; % mean Calpha-Calpha distance for direct neighbor residues


switch ENM_param.parametrization
    case 'cutoff10'
        ENM_param.p_ANM=0; %  exponent for distance-dependent force constants ANM
        ENM_param.cutoff_const_ANM=10;
        ENM_param.cutoff_log_ANM=0;
        ENM_param.connect=0;
        ENM_param.C_connect=0;
        ENM_param.p_connect=0;
        ENM_param.C_cart=6;
    case 'cutoff13'
        ENM_param.p_ANM=0; %  exponent for distance-dependent force constants ANM
        ENM_param.cutoff_const_ANM=13;
        ENM_param.cutoff_log_ANM=0;
        ENM_param.connect=0;
        ENM_param.C_connect=0;
        ENM_param.p_connect=0;
        ENM_param.C_cart=6;
    case 'Hinsen'
        ENM_param.p_ANM=6; %  exponent for distance-dependent force constants ANM
        ENM_param.cutoff_const_ANM=1e6;
        ENM_param.cutoff_log_ANM=0;
        ENM_param.connect=1;
        ENM_param.C_connect=(8.6e5*nm2A^(-3)*rCaCa-2.39e5*nm2A^(-2))/RT;
        ENM_param.p_connect=100;
        ENM_param.C_cart=(128*nm2A^(4))/RT;
    case 'Jeschke'
        ENM_param.p_ANM=6; %  exponent for distance-dependent force constants ANM
        ENM_param.cutoff_const_ANM=1e6;
        ENM_param.cutoff_log_ANM=0;
        ENM_param.connect=2;
        ENM_param.C_connect=19.62*kcal2kJ;
        ENM_param.p_connect=1.789;
        ENM_param.C_cart=6*kcal2kJ;
    case 'ed-ENM'
        ENM_param.p_ANM=6; %  exponent for distance-dependent force constants ANM
        ENM_param.cutoff_const_ANM=-2.9;
        ENM_param.cutoff_log_ANM=2.9;
        ENM_param.connect=3;
        ENM_param.C_connect=(2.5e4*nm2A^(-2))/RT;
        ENM_param.p_connect=2;
        ENM_param.C_cart=(19.534*nm2A^(4))/RT;
    case 'ed-ENM-p'
        ENM_param.p_ANM=6; %  exponent for distance-dependent force constants ANM
        ENM_param.cutoff_const_ANM=1e6;
        ENM_param.cutoff_log_ANM=0;
        ENM_param.connect=3;
        ENM_param.C_connect=60*kcal2kJ;
        ENM_param.p_connect=2;
        ENM_param.C_cart=6*kcal2kJ;
end

