function [C,R0,nu,raxis,Rg_distr,P,R0_P,nu_P,all_k,all_Rg,all_R2] = local_compactness(coor,N,n,p)
% [C,R0,nu,raxis,Rg_distr,P,R0_P,nu_P] = local_compactness(coor,N,n)
% 
% Local compactness analysis for an ensemble of N Calpha coils with n
% residues each
%
% coor      (N*n,3) Cartesian coordinate array
% N         number of structures in ensemble
% n         number of residues (Calpha atoms) in one structure
% p         population vector (1,N), if empty or missing, all populations
%           are 1/N, otherwise populations are normalized
%
% C         local compactness matrix
% R0        segment length in Rg = R0 n`nu
% nu        scaling (chain extension) exponent in Rg = R0 n^nu
% raxis     distance axis for Rg distribution matrix
% Rg_distr  matrix of radius of gyration distribution versus segment length
% P         proximity matrix
% R0_P      sement length in sqrt(<R^2>) = R0_P n^nu_P
% nu_P      scaling (chain extension) exponent in sqrt(<R^2>) = R0_P n^nu_P
%
% G. Jeschke, 2019-2020

if ~exist('p','var') || isempty(p)
    p = ones(1,N)/N;
else
    p = p/sum(p);
end

if ~exist('options','var') || isempty(options) || ~isfield(options,'verbose')
    options.verbose = false;
end

if ~isfield(options,'chain') || isempty(options.chain)
    options.chain = 'A';
end

R0 = 2; % starting value for effective segment length
nu = 0.6; % starting value for Flory exponent
increments_per_Angstroem = 10;

[m,c] = size(coor);
if m ~= N*n
    warning('local_compactness: coordinate matrix size does not match ensemble size and residue number');
    return
end
if c ~= 3
    warning('local_compactness: second dimension of coordinate array should be 3');
end

C = zeros(n);
P = zeros(n);
fully_extended = zeros(n,3);
fully_extended(:,1) = linspace(0,(n-1)*3.8,n)';
rmax = radius_of_gyration(fully_extended); % radius of gyration of the completely extended chain
raxis = linspace(0,ceil(rmax),1+(increments_per_Angstroem*ceil(rmax)));
Rg_distr = zeros(length(raxis),n);
Rg_distr(1,1) = 1; % chain with length zero has radius of gyration of zero

exp_rmax = 0; % initailize experimentally observed maximum radius of gyration
all_Rg = zeros(1,n*(n-1)/2);
all_R2 = zeros(1,n*(n-1)/2);
all_Rg_poi = 0;
all_k = all_Rg;

% double loop over all residue pairs, indices as in the white paper on
% ensemble analysis
for i = 1:n-1
    for j = i+1:n
        all_Rg_poi = all_Rg_poi + 1;
        % extract coordinates for this segment from all structures
        l = j - i + 1;
        k = l - 1;
        segment_coor = zeros(N*l,3);
        seg_ind = 1:l;
        orig_ind = i:j;
        R2 = 0;
        for s = 1:N
            R2 = R2 + p(s)*norm(coor(orig_ind(end),:)-coor(orig_ind(1),:))^2;
            segment_coor(seg_ind,:) = coor(orig_ind,:);
            seg_ind = seg_ind + l;
            orig_ind = orig_ind + n;
        end
        P(i,j) = sqrt(R2);
        P(j,i) = sqrt(R2);
        all_R2(all_Rg_poi) = sqrt(R2);
        Rg = ensemble_radius_of_gyration(segment_coor,N,p);
        C(i,j) = Rg;
        C(j,i) = Rg;
        all_Rg(all_Rg_poi) = Rg;
        all_k(all_Rg_poi) = k;
        rpoi = 1 + round(increments_per_Angstroem * Rg);
        Rg_distr(rpoi,k) = Rg_distr(rpoi,k) + 1;
        if Rg > exp_rmax
            exp_rmax = Rg;
        end
    end
end
% cut histograms to maximal encountered radius of gyration 
maxrpoi = ceil(increments_per_Angstroem * exp_rmax);
Rg_distr = Rg_distr(1:maxrpoi,:);
raxis = raxis(1:maxrpoi);
% normalize histograms
for k = 1:n
    Rg_distr(:,k) = Rg_distr(:,k)/sum(Rg_distr(:,k));
end
% fit random coil parameters
v0(1) = R0;
v0(2) = nu;
v = fminsearch(@msd_flory,v0,[],all_k,all_Rg);
R0 = v(1);
nu = v(2);

% fit random coil parameters for mean square end-to-end distance
v0(1) = sqrt(6)*R0;
v0(2) = nu;
v = fminsearch(@msd_flory,v0,[],all_k,all_R2);
R0_P = v(1);
nu_P = v(2);

% Subtract reference state from compactness matrix and normalize to
% reference state
% double loop over all residue pairs, indices as in the white paper on
% ensemble analysis
for i = 1:n-1
    for j = i+1:n
        k = j-i;
        Rg = R0*k^nu;
        C(i,j) = (C(i,j)-Rg)/Rg;
        C(j,i) = C(i,j);
        R2 = R0_P*k^nu_P;
        P(i,j) = (P(i,j)-R2)/R2;
        P(j,i) = P(i,j);
    end
end

function msd = msd_flory(v,kvec,Rgvec)

sim_Rg = v(1)*kvec.^v(2);
msd = sum(Rgvec-sim_Rg).^2/length(Rgvec);

