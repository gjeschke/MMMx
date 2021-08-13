function [Hessian,cutoff] = get_ANM_Hessian(CA_model,ENM_param)
% [Hessian,cutoff] = GET_ANM_HESSIAN(CA_model,ENM_param)
%
% returns the Hessian for a fully connected anisotropic network model (ANM)
% when given the Calpha coordinates
% the force constants have inverse 6th power dependence on interresidue
% distance, except for next neighbors (i,i+1) and second neighbors (i,i+2)
% along the polypetide chain
% parameters of the elastic network model are loaded from 
%
% see: I. Bahar, T. R. Lezon, A. Bakan, I. H. Shrivastava,
% Normal Mode Analysis of Biomolecular Structures: Functional Mechanisms of 
% Membrane Proteins. 
% Chem. Rev. 110 (2010) 1463-1497.
% section 2.3.1
% and:
% L. Yang, G. Song, R L. Jernigan, Proc. Natl. Acad. Sci. 106 (2009)
% 12347-12352.
% and:
% L. Orellana, M. Rueda, C. Ferrer-Costa, J. R. Lopez-Blanco, P. Chacon, M.
% Orozco, J. Chem. Theory Comput. 2010, 6, 2910-2923
%
% Input:
%
% CA_model  struct with information on the Calpha trace model, fields
%           .coor       array(N,3) of Cartesian coordinates of Calpha atoms
%           .masses     array(N,1) of residue masses
%           .bfactors   array(N,1) with isotropic Bfactors, NaN if missing
%           .chains     array(N,1) of chain indices
%           .addresses  cell(N,1) of corresponding residue addresses
% ENM_param elastic network parameters, default to output of get_ENM.m
%
%
% parameters are defined in
% initialize_MMM.m  and stored in global variable ENM_param:
%      .p_ANM               exponent for distance-dependent force constant
%      .cutoff_const_ANM    constant contribution to cutoff distance
%      .cutoff_log_ANM      logarithmic scaling of cutoff with number of
%                           residues
%      .connect             number of connected neighbors with larger force
%                           constants
%      .C_connect           connected force constant for direct neighbor in
%                           kcal/?^2
%      .p_connect           exponent for inverse scaling of connected force
%                           constant with distance in residue sequence
%      .C_cart              unconnected force constant for direct neighbor 
%                           in kcal/?^2
%
% an empty matrix is returned if the routine fails (no connection between
% residues)
%
% G. Jeschke, 2010-2021

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% G. Jeschke, 2010-2021

% define the elastic network model
if ~exist('ENM_param','var') || isempty(ENM_param)
    ENM_param = get_ENM;
end

min_cutoff=7.0; % minimum cutoff distance

[m,~] = size(CA_model.coor); % m is number of Calpha atoms

% cutoff for Cartesian distances
cutoff = ENM_param.cutoff_log_ANM*log(m) + ENM_param.cutoff_const_ANM;
if cutoff < min_cutoff
    cutoff = min_cutoff;
end

% make vector of bond forces
bond_forces = zeros(1,ENM_param.connect);
for k = 1:ENM_param.connect
    bond_forces(k) = ENM_param.C_connect/(k^ENM_param.p_connect);
end

% make membrane scaling matrix for imANM by Lezon/Bahar
%
if ENM_param.imANM
    sqrt_gz = sqrt(ENM_param.imANM_scaling);
    anisotropy = ones(1,1);
    anisotropy(1,3) = sqrt_gz;
    anisotropy(2,3) = sqrt_gz;
    anisotropy(3,1) = sqrt_gz;
    anisotropy(3,2) = sqrt_gz;
    anisotropy(3,3) = ENM_param.imANM_scaling;
end

Hessian=zeros(3*m);
for i=1:m-1
    basi=3*(i-1);
    for j=i+1:m
        basj=3*(j-1);
        rvec = CA_model.coor(i,:) - CA_model.coor(j,:); % distance vector of Calpha atoms
        r=norm(rvec);
        if r<=cutoff
            if ENM_param.mass_weighting
                m1 = CA_model.masses(i);
                m2 = CA_model.masses(j);
                wgamma = 1/sqrt(m1*m2);
            else
                wgamma = 1;
            end
            if CA_model.chains(i) == CA_model.chains(j) == 0 % both C_alpha in the same chain
                num1 = get_residue_number(CA_model.addresses{i});
                num2 = get_residue_number(CA_model.addresses{j});
                diff = num2-num1;
                if diff <= ENM_param.connect && r <= ENM_param.bond_max(diff)
                    cgamma = -bond_forces(diff);
                else
                    cgamma = -ENM_param.C_cart * r^(-ENM_param.p_ANM);
                end
            else
                cgamma = -ENM_param.C_cart * r^(-ENM_param.p_ANM);
            end
            submat = kron(rvec,rvec')/r^2; % see Eq. (24) I. Bahar et al. Chem. Rev.
            if ENM_param.imANM
                submat=submat.*anisotropy;
            end
            Hessian(basi+1:basi+3,basj+1:basj+3) = wgamma*cgamma*submat;
            Hessian(basj+1:basj+3,basi+1:basi+3) = wgamma*cgamma*submat;
        end
    end
end

conn=0;
for i = 1:m
    basi = 3*(i-1);
    submat=zeros(3);
    for j = 1:m
        basj = 3*(j-1);
        submat = submat + Hessian(basi+1:basi+3,basj+1:basj+3);
    end
    conn = conn + trace(submat);
    Hessian(basi+1:basi+3,basi+1:basi+3) = -submat;
end

if conn == 0 % ANM is unconnected
    Hessian=[];
end

function num = get_residue_number(address)

poi = strfind(address,')');
if isempty(poi)
    poi = 0;
end
num = str2double(address(poi+1:end));
