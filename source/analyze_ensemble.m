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
%               .mono (1,N) double, numbers of the residues
%               .bb   (1,C) cell of (N,3) coordinate array for backbone
%                     atoms
%               .slc  single-letter code for residues
% pop           (1,C) population vector for C conformers, defaults to
%               uniform populations
% options       options for output to be computed, struct with fields
%               .chain_mode     Boolean, if true, analysis for individual
%                               chains, defaults to false
%               .Rg             Boolean, if true, radius of gyration is
%                               analyzed, defaults to true
%               .pair_rmsd      Boolean, pairwise root mean square 
%                               deviation,
%                               defaults to false
%               .pair_drms      Boolean, pairwise distance root mean square 
%                               deviation, ensemble width, and density,
%                               defaults to true
%               .superimpose    Boolean, optimal superposition of
%                               conformers for computing pair rmsd,
%                               defaults to true, effect only if .pair_rmsd
%                               is true
%               .sorted         Boolean, recursive hierarchical clustering
%                               sorting, defaults to false, forces
%                               .pair_rmsd to true if neither .pair_rmsd
%                               nor .pair.drms is set, if .pair_drms is
%                               set, sorting is by distance root mean
%                               square deviation
%               .pair_corr      Boolean, if true, compute residue pair
%                               correlations, defaults to true 
%               .compactness    Boolean, if true, compactness analysis is
%                               performed, defaults to true
%               .type_analysis  Boolean, if true, proximity analysis is
%                               performed per residue type, requires
%                               .compactness to be true and chain must be a
%                               biopolymer, defaults to true if .chain_mode
%                               is true
%
% OUTPUT
% measures      ensemble measures, empty if no input information
%               struct with fields for all chains (chain mode) or a single
%               field .all (chain mode false), each field has subfields
%               .Rg         radius of gyration
%               .Rg_std     standard deviation of radius of gyration
%               .Rg_all     radii of gyration for all conformers
%               .pair_rmsd  matrix (C,C) of pair rmsd between conformers
%               .width      ensemble width (population-weighted)
%               .density    mean distance (rmsd) between next-neighbor
%                           conformers 
%               .sorting    indices of the sorted ensemble
%               .R0_g       segment length for Flory fit of radius of
%                           gyration
%               .nu_g       scaling exponent for Flory fit of radius of
%                           gyration
%               .Rg_axis    distance axis for radius of gyration
%                           distribution
%               .Rg_distr   radius of gyration distribution
%               .R0_ee      segment length for Flory fit of mean square
%                           end-to-end distance
%               .nu_ee      scaling exponent for Flory fit of mean square
%                           end-to-end distance
%               .seg_len    vector of the lengths of all segments
%               .all_Rg     vector of the radii of gyration of all segments
%               .all_R2     vector of mean square end-to-end distances of
%                           all segments
%               depending on options, not all fields may be present
% correlations  correlation matrices, empty if no input information
%               struct with fields for all chains (chain mode) or a single
%               field .all (chain mode false), each field has subfields
%               .pair_axis  cell string, axis ('chain.residue') for pair
%                           correlations
%               .pair_corr  double (N,N) pair correlation matrix 
%                           (standard deviation of pair distance)
%               .dmat       double (N,N) mean distance matrix
%               .compact    double(N,N) compactness matrix
%               .proximity  double(N,N) proximity matrix
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

% set default options
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

if ~isfield(options,'pair_rmsd') || isempty(options.pair_rmsd)
    options.pair_rmsd = false;
end

if ~isfield(options,'superimpose') || isempty(options.superimpose)
    options.superimpose = true;
end

if ~isfield(options,'pair_drms') || isempty(options.pair_drms)
    options.pair_drms = true;
end

if ~isfield(options,'sorted') || isempty(options.sorted)
    options.sorted = false;
end

if options.sorted && ~options.pair_drms
    options.pair_rmsd = true;
end

if ~isfield(options,'pair_corr') || isempty(options.pair_corr)
    options.pair_corr = true;
end

if ~isfield(options,'compactness') || isempty(options.compactness)
    options.compactness = true;
end

if ~isfield(options,'type_analysis') || isempty(options.type_analysis)
    if options.chain_mode
        options.type_analysis = true;
    else
        options.type_analysis = false;
    end
end

if options.type_analysis
    options.chain_mode = true;
end

if ~options.chain_mode
    % determine total size of coordinate array
    N = 0;
    for kc = 1:length(chains)
        bb0 = backbones.(chains{kc}).bb{1};
        [N0,~] = size(bb0);
        N = N + N0;
    end
    backbones.all.chains = blanks(N);
    backbones.all.mono = zeros(1,N);
    backbones.all.slc = blanks(N);
    poi = 0;
    for kc = 1:length(chains)
        for km = 1:length(backbones.(chains{kc}).mono)
            backbones.all.chains(poi+km) = chains{kc};
            backbones.all.slc(poi+km) = backbones.(chains{kc}).slc(km);
            backbones.all.mono(poi+km) = backbones.(chains{kc}).mono(km);
        end
        poi = poi + length(backbones.(chains{kc}).mono);
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
    chains = fieldnames(backbones); % assign 'all' to chain list
end

for kc = 1:length(chains)
    % Radius of gyration
    if options.Rg
        Rg_all = zeros(size(pop));
        for c = 1:C
            coor = backbones.(chains{kc}).bb{c};
            Rg_all(c) = gyration_radius(coor);
        end
        measures.(chains{kc}).Rg = sqrt(sum(pop.*Rg_all.^2)); % we average the square of the radius of gyration
        measures.(chains{kc}).Rg_all = Rg_all;
        measures.(chains{kc}).Rg_std = std(Rg_all);
    end
    % pair rmsd, ensemble width, and ensemble density
    if options.pair_rmsd
        pair_rmsd = zeros(C);
        for c1 = 1:C-1
            coor1 = backbones.(chains{kc}).bb{c1};
            for c2 = c1+1:C
                coor2 = backbones.(chains{kc}).bb{c2};
                if options.superimpose
                    rmsd = rmsd_superimpose(coor1,coor2);
                else
                    [A,~] = size(coor1);
                    rmsd = sqrt(sum(sum((coor1-coor2).^2))/A);
                end
                pair_rmsd(c1,c2) = rmsd;
                pair_rmsd(c2,c1) = rmsd;
            end
        end
        measures.(chains{kc}).pair_rmsd = pair_rmsd;        
        if options.sorted && ~options.pair_drms
            rmsdsum = sum(pair_rmsd);
            [~,sorting] = sort(rmsdsum);
            pair_rmsd = pair_rmsd(sorting,sorting);
            [new_pair_rmsd,new_ordering] = cluster_sorting(pair_rmsd,pop,1,sorting);
            measures.(chains{kc}).sorting = new_ordering;
            measures.(chains{kc}).pair_rmsd = new_pair_rmsd;
        else
            measures.(chains{kc}).sorting = 1:C;
        end
    end
    if options.pair_drms
        pair_drms = zeros(C);
        coor0 = backbones.(chains{kc}).bb{1};
        [N_atoms,~] = size(coor0);
        N2 = N_atoms*(N_atoms-1)/2;
        for c1 = 1:C-1
            coor1 = backbones.(chains{kc}).bb{c1};
            for c2 = c1+1:C
                coor2 = backbones.(chains{kc}).bb{c2};
                dmat1 = coor2dmat(coor1); % distance matrix for first conformer
                dmat2 = coor2dmat(coor2); % distance matrix for second conformer
                drms = sum(sum((dmat1-dmat2).^2)); % mean-square deviation for this conformer pair
                pair_drms(c1,c2) = drms;
                pair_drms(c2,c1) = drms;
            end
        end
        pair_drms = sqrt(pair_drms/N2); % normalized root mean square deviation
        measures.(chains{kc}).pair_drms = pair_drms;
        measures.(chains{kc}).width = get_ensemble_width(pair_drms,pop);
        measures.(chains{kc}).density = get_ensemble_density(pair_drms);
        if options.sorted
            rmsdsum = sum(pair_drms);
            [~,sorting] = sort(rmsdsum);
            pair_drms = pair_drms(sorting,sorting);
            [new_pair_drms,new_ordering] = cluster_sorting(pair_drms,pop,1,sorting);
            measures.(chains{kc}).sorting = new_ordering;
            measures.(chains{kc}).pair_drms = new_pair_drms;
        else
            measures.(chains{kc}).sorting = 1:C;
        end
    end
    % pair correlations
    if options.pair_corr
        [corr_isotropic,corr_axis,dmat] = get_correlations(backbones,pop);
        correlations.(chains{kc}).pair_corr = corr_isotropic;
        correlations.(chains{kc}).pair_axis = corr_axis;
        correlations.(chains{kc}).dmat = dmat;
    end
    % compactness analysis
    if options.compactness
        if options.type_analysis
            % load monomer definitions
            defs = load('monomers.mat');
            if backbones.(chains{kc}).type == 1 % peptide
                options.slc = defs.residue_slc(1:20);
            elseif backbones.(chains{kc}).type == 2 % nucleic acid
                options.slc = defs.residue_slc(21:28);
            else % should not happen, but if not a biopolymer, switch this analysis off
                options.type_analysis = false;
            end
            [Comp,R0,nu,raxis,Rg_distr,Prox,R0_P,nu_P,all_k,all_Rg,all_R2,mean_R2,R0_seglen,nu_seglen,Rdev,type_proxi]...
                = local_compactness(backbones.(chains{kc}),pop,options);
            correlations.(chains{kc}).type_proxi = type_proxi;
            correlations.(chains{kc}).res_types = options.slc;
        else % no type analysis
            [Comp,R0,nu,raxis,Rg_distr,Prox,R0_P,nu_P,all_k,all_Rg,all_R2,mean_R2,R0_seglen,nu_seglen,Rdev]...
                = local_compactness(backbones.(chains{kc}),pop,options);
        end
        % make a more compact version of sequence length distribution
        seg_lengths = 1:max(all_k);
        min_Rg = 1e6*ones(size(seg_lengths));
        max_Rg = zeros(size(seg_lengths));
        min_R2 = 1e6*ones(size(seg_lengths));
        max_R2 = zeros(size(seg_lengths));
        for j = 1:length(all_k)
            k = all_k(j);
            if all_Rg(j) < min_Rg(k)
                min_Rg(k) = all_Rg(j);
            end
            if all_Rg(j) > max_Rg(k)
                max_Rg(k) = all_Rg(j);
            end
            if all_R2(j) < min_R2(k)
                min_R2(k) = all_R2(j);
            end
            if all_R2(j) > max_R2(k)
                max_R2(k) = all_R2(j);
            end
        end
        measures.(chains{kc}).R0_g = R0;
        measures.(chains{kc}).nu_g = nu;
        measures.(chains{kc}).Rg_axis = raxis;
        measures.(chains{kc}).Rg_distr = Rg_distr;
        measures.(chains{kc}).R0_ee = R0_P;
        measures.(chains{kc}).nu_ee = nu_P;
        measures.(chains{kc}).seg_len = all_k;
        measures.(chains{kc}).all_Rg = all_Rg;
        measures.(chains{kc}).all_R2 = all_R2;
        measures.(chains{kc}).mean_R2 = mean_R2;
        measures.(chains{kc}).R0_seglen = R0_seglen;
        measures.(chains{kc}).nu_seglen = nu_seglen;
        measures.(chains{kc}).seg_lengths = seg_lengths;
        measures.(chains{kc}).min_Rg = min_Rg;
        measures.(chains{kc}).max_Rg = max_Rg;
        measures.(chains{kc}).min_R2 = min_R2;
        measures.(chains{kc}).max_R2 = max_R2;
        correlations.(chains{kc}).compact = Comp;
        correlations.(chains{kc}).proximity = Prox;        
        correlations.(chains{kc}).Rdev = Rdev;        
    end
end

function sigma = get_ensemble_width(pair_rmsd,pop)

[C,~] = size(pair_rmsd);

msdsum = 0;
popsum = 0;
for c1 = 1:C-1
    for c2 = c1+1:C
        msdsum = msdsum + pop(c1)*pop(c2)*pair_rmsd(c1,c2)^2;
        popsum = popsum + pop(c1)*pop(c2);
    end
end
sigma = sqrt(msdsum/popsum);

function density = get_ensemble_density(pair_rmsd)

[C,~] = size(pair_rmsd);
msq = zeros(1,C);
for c = 1:C
    cmsd = pair_rmsd(c,:).^2;
    cmsd(c) = 1e12;
    msq(c) = min(cmsd);
end
density = sqrt(mean(msq));

function [corr_isotropic,corr_axis,dmat] = get_correlations(backbones,pop)

chains = fieldnames(backbones);
nch = length(chains);
C = length(backbones.(chains{1}).bb);
resnums = zeros(1,nch);
for kc = 1:nch
    bb0 = backbones.(chains{kc}).bb{1};
    [nr,~] = size(bb0);
    resnums(kc) = nr;
end
resnum = sum(resnums);
corr_isotropic = zeros(resnum);
dmat = zeros(resnum);
corr_axis = cell(1,resnum);
poi1 = 0;
dmats = cell(1,C);
for kc = 1:nch
    chain = chains{kc};
    for kr = 1:resnums(kc)
        poi1 = poi1 + 1;
        if strcmp(chain,'all')
            corr_axis{poi1} = sprintf('%s.%i',backbones.all.chains(kr),backbones.all.mono(kr));
        else
            corr_axis{poi1} = sprintf('%s.%i',chain,backbones.(chain).mono(kr));
        end
    end
end
for c = 1:C
    coor = zeros(resnum,3);
    poi2 = 0;
    for kc = 1:nch
        coor(1+poi2:resnums(kc)+poi2,:) = backbones.(chain).bb{c};
        poi2 = poi2 + resnums(kc);
    end
    a=coor';
    ab = coor*a;
    aa=sum(a.*a,1); 
    n =size(aa,2);
    aa2 = repmat(aa,[n 1]);
    dmats{c} = sqrt(abs(aa2' + aa2 - 2*ab));
end
for r1 = 1:resnum
    for r2 = 1:resnum
        rvec = zeros(size(pop));
        for c = 1:C
            rvec(c) = dmats{c}(r1,r2);
        end
        rmean = mean(rvec);
        rvec = rvec - rmean;
        corr_isotropic(r1,r2) = sqrt(sum(pop.*rvec.^2));
        corr_isotropic(r2,r1) = corr_isotropic(r1,r2);
        dmat(r1,r2) = rmean;
        dmat(r2,r1) = rmean;
    end
end

function [Comp,R0,nu,raxis,Rg_distr,Prox,R0_P,nu_P,all_k,all_Rg,all_R2,mean_R2,R0_seglen,nu_seglen,Rdev,type_proxi] = local_compactness(backbone,pop,options)
% [Comp,R0,nu,raxis,Rg_distr,P,R0_P,nu_P] = local_compactness(backbones,pop)
% 
% Local compactness analysis for an ensemble of N Calpha coils with n
% residues each
%
% Input:
%
% backbone      backbone structure fir one chain with fields
%               .type 1 for peptide, 2 for nucleic acid, for mixed chains,
%                     the type with more backbone atoms prevails
%               .mono (1,nm) double, numbers of the residues
%               .bb   (1,C) cell of (nm,3) coordinate array for backbone
%                     atoms
%               .slc  single-letter code for residues
% pop           (1,C) population vector for C conformers, defaults to
%               uniform populations
% options       option struct with fields
%               .type_analysis  Boolean, if true, proximity analysis is
%                               performed per residue type, defaults to
%                               false
%               .slc            single_letter codes for indexing
%
% Output:
%
% Comp      local compactness matrix
% R0        segment length in Rg = R0 n`nu
% nu        scaling (chain extension) exponent in Rg = R0 n^nu
% raxis     distance axis for Rg distribution matrix
% Rg_distr  matrix of radius of gyration distribution versus segment length
% Prox      proximity matrix
% R0_P      segment length in sqrt(<R^2>) = R0_P n^nu_P
% nu_P      scaling (chain extension) exponent in sqrt(<R^2>) = R0_P n^nu_P
% all_k     vector of the lengths of all encountered segments
% all_Rg    vector of the radii of gyration of all encountered segments
% all_R2    vector of mean square end-to-end distances of all encountered
%           segments 
% mean_R2   k-averaged end-to-end distance
% R0_seglen segment length in fit to mean_R2
% nu_seglen scaling exponent in fit to mean R2
% Rdev      matrix with deviation of end-to-end distances from mean_R2
%
% G. Jeschke, 2019-2020

% empty type analysis output
type_proxi = [];

R0 = 2; % starting value for effective segment length
nu = 0.6; % starting value for Flory exponent
increments_per_Angstroem = 10; % resolution

C = length(backbone.bb);

if ~exist('pop','var') || isempty(pop)
    pop = ones(1,C)/C;
else
    pop = pop/sum(pop);
end

if ~exist('options','var') || isempty(options) || ~isfield(options,'type_analysis')
    options.type_analysis = false;
end

coor = backbone.bb{1};
[N,~] = size(coor); % number of residues

Comp = zeros(N); % initalize compactness matrix
Prox = zeros(N); % initialize proximity matrix
fully_extended = zeros(N,3);
fully_extended(:,1) = linspace(0,(N-1)*3.8,N)';
rmax = gyration_radius(fully_extended); % radius of gyration of the completely extended chain
raxis = linspace(0,ceil(rmax),1+(increments_per_Angstroem*ceil(rmax))); % define distance axis
Rg_distr = zeros(length(raxis),N); % radii of gyration distribution per segment length
Rg_distr(1,1) = 0; % chain with length zero has radius of gyration of zero

exp_rmax = 0; % initialize experimentally observed maximum radius of gyration
all_Rg = zeros(1,N*(N-1)/2); % all encountered radii of gyration
all_R2 = zeros(1,N*(N-1)/2); % all encountered mean square end-to-end distances
all_Rg_poi = 0; % pointer for radii of gyration
all_k = all_Rg; % sgement length
all_i = all_Rg; % index of first residue
all_j = all_Rg; % index of second residue
if options.type_analysis
    all_types = zeros(N*(N-1)/2,2); % all residue-type pairs
end

% double loop over all residue pairs, indices as in the white paper on
% ensemble analysis
for i = 1:N-1 % first residue of segment
    for j = i+1:N % last residue of segment
        all_Rg_poi = all_Rg_poi + 1; % increment counter
        if options.type_analysis % assign residue-type pair, if necessary
            type_i = strfind(options.slc,backbone.slc(i));
            if ~isempty(type_i) % only standard residues are typed
                all_types(all_Rg_poi,1) = type_i;
            end
            type_j = strfind(options.slc,backbone.slc(j));
            if ~isempty(type_j)
                all_types(all_Rg_poi,2) = type_j; 
            end
        end
        % extract coordinates for this segment from all structures
        len = j - i + 1; % length of segment plus one (Zaunlattenprinzip)
        k = len - 1;
        segment_coor = zeros(C*len,3); % segment coordinates for all conformers
        seg_ind = 1:len;  % indices of residues in segment
        orig_ind = i:j; % indices of segment in complete chain
        R2 = 0;
        for c = 1:C
            coor = backbone.bb{c};
            R2 = R2 + pop(c)*norm(coor(orig_ind(end),:)-coor(orig_ind(1),:))^2; % means square end-to-end distance of segment
            segment_coor(seg_ind,:) = coor(orig_ind,:); % coordinates for whole segment
            seg_ind = seg_ind + len; % increment coordinate index by segment length
        end
        % Proximity matrix element
        Prox(i,j) = sqrt(R2);
        Prox(j,i) = sqrt(R2);
        all_R2(all_Rg_poi) = sqrt(R2); % vector of all root mean square end-to-end distances
        % compute ensemble radius of gyration
        Rg2 = 0;
        for c = 1:C
            coor = segment_coor(1+(c-1)*len:c*len,:); % segment coordinates for conformer c
            rc = mean(coor); % center coordinate
            rel_coor = coor - repmat(rc,len,1); % subtract center coordinate
            Rg2 = Rg2 + pop(c)*sum(sum(rel_coor.^2))/len; % mean square radiues of gyration
        end
        Rg = sqrt(Rg2); % root mean square radius of gyration
        % Compactness matrix element
        Comp(i,j) = Rg;
        Comp(j,i) = Rg;
        all_Rg(all_Rg_poi) = Rg; % vector of all radii of gyration
        all_k(all_Rg_poi) = k; % vector of all segment lengths
        all_i(all_Rg_poi) = i; % vector of all indices of residue 1
        all_j(all_Rg_poi) = j; % vector of all indices of residue 2
        rpoi = 1 + round(increments_per_Angstroem * Rg); % pointer into radii of gyration histogram
        Rg_distr(rpoi,k) = Rg_distr(rpoi,k) + 1; % update histogram for this segment length
        if Rg > exp_rmax % update maximum encountered radius of gyration
            exp_rmax = Rg;
        end
    end
end
% cut histograms and distance axis to maximal encountered radius of gyration 
maxrpoi = ceil(increments_per_Angstroem * exp_rmax);
Rg_distr = Rg_distr(1:maxrpoi,:);
raxis = raxis(1:maxrpoi);
% normalize histograms
for k = 1:N
    Rg_distr(:,k) = Rg_distr(:,k)/sum(Rg_distr(:,k));
end
% fit random coil parameters
v0(1) = R0;
v0(2) = nu;
v = fminsearch(@msd_flory,v0,[],all_k,all_Rg); % not the best minimizer, but the problem is simple
R0 = v(1);
nu = v(2);

% fit random coil parameters for mean square end-to-end distance
v0(1) = sqrt(6)*R0;
v0(2) = nu;
v = fminsearch(@msd_flory,v0,[],all_k,all_R2); % not the best minimizer, but the problem is simple
R0_P = v(1);
nu_P = v(2);


% Subtract reference state from compactness and prioximity matrix and normalize to
% reference state
% double loop over all residue pairs, indices as in the white paper on
% ensemble analysis
for i = 1:N-1
    for j = i+1:N
        k = j-i;
        Rg = R0*k^nu;
        Comp(i,j) = (Comp(i,j)-Rg)/Rg;
        Comp(j,i) = Comp(i,j);
        R2 = R0_P*k^nu_P;
        Prox(i,j) = (Prox(i,j)-R2)/R2;
        Prox(j,i) = Prox(i,j);
    end
end

kaxis = 1:max(all_k);
mean_R2 = zeros(size(kaxis));
for k = kaxis
    these_R2 = all_R2(all_k == k);
    mean_R2(k) = mean(these_R2);
end

% fit random coil parameters to the mean curve
v0(1) = sqrt(6)*R0;
v0(2) = nu;
v = fminsearch(@msd_flory,v0,[],kaxis,mean_R2); % not the best minimizer, but the problem is simple
R0_seglen = v(1);
nu_seglen = v(2);



Rdev = zeros(size(Prox));
for me = 1:length(all_R2) % loop over all elements of the upper triangular matrix
    Rdev(all_i(me),all_j(me)) = (all_R2(me) - mean_R2(all_k(me))); %/all_R2(me); 
    Rdev(all_j(me),all_i(me)) = Rdev(all_i(me),all_j(me));
end


% do type analysis if requested
if options.type_analysis % assign residue-type pair, if necessary
    type_proxi = zeros(length(options.slc)); 
    type_count = zeros(length(options.slc)); 
    for pair = 1:N*(N-1)/2 % loop over all residue pairs
        type_i = all_types(pair,1);
        type_j = all_types(pair,2);
        if type_i ~= 0 && type_j ~= 0 % only pairs with two native residues
            R2 = R0_P*all_k(pair)^nu_P; % expected distance from compactness
            type_count(type_i,type_j) = type_count(type_i,type_j) + 1;
            type_count(type_j,type_i) = type_count(type_j,type_i) + 1;
            type_proxi(type_i,type_j) = type_proxi(type_i,type_j) + all_R2(pair)/R2;
            type_proxi(type_j,type_i) = type_proxi(type_j,type_i) + all_R2(pair)/R2;
        end
    end
    type_proxi = type_proxi./type_count;
end

function msd = msd_flory(v,kvec,Rgvec)

sim_Rg = v(1)*kvec.^v(2);
msd = sum(Rgvec-sim_Rg).^2/length(Rgvec);

