function [fom,rmsd,network,dxmat,ddr,diagnostics] = fit_by_ANM_thermal(network0,ddr,options,CA_model,ENM_param)
% [fom,rmsd,network,dxmat,ddr,best,diagnostics] = FIT_BY_ANM_THERMAL(network0,ddr,options,CA_model,ENM_param)
%
% fits a network model to distance distribution restraints based on 
% equipartitioning of energy over the modes, mode signs are derived
% from direction of the improvement in distance restraint r.m.s.d.
%
% specifics of local structure restraints to fix C_alpha-C_alpha distances
% or to stabilize secondary structure are defined in ENM_param.fix_local
% see get_ENM.m for initialization of ENM_param
%
% Input:
% network0  matrix [n,3] of coordinates of network points (knots),
%           restraint variable ddr must refer labelled sites to network
%           knots (C_alpha coordinates of residues), these references are
%           given for the k-th restraint by ddr(k).index1 and ddr(k).index2
%           indexing of residues in network is supposed to conform to
%           CA_model
% ddr       distance distribution restraints, (1,n) struct for n restraints
%           all in  units of Angstroem
%           .r          mean distances
%           .sigr       uncertainties of mean distances
%           .xyz1       Cartesian coordinates of mean position first label
%           .xyz2       Cartesian coordinates of mean position second label
%           .index1     number of first residue in Calpha network model
%           .index2     number of second residue in Calpha network model
% options   structure with fit options (optional, fields also optional)
%           .maxbas     maximum basis size for adaptive basis extension,
%                       if empty or missing, basis extends to 1/3 of all
%                       modes within the maximum number of iteration cycles
%                       ENM_param.cycles
%           .sat_cyc    iteration cycle, at which maximum basis size is 
%                       attained, defaults to 100, is unused if .maxbas is
%                       missing or empty
% CA_model  Calpha model for m residues, struct with fields
%           .coor       (m,3) double Cartesian coordinates of CA atoms
%           .masses     (m,1) double, residue masses in g/mol (unused)
%           .bfactors   (m,1) double, B factors (unused, can be NaN)
%           .chains     (m,1) double, chain indices
%           .addresses  (m,1) cell of strings with residue addresses
% ENM_param parameters of the elastic network model, see get_ENM.m
%
% Output:
% fom           figure of merit, quantifies deviation of the fitted model
%               from restraints, maximum likelihood estimated based on mean
%               distance and standard deviation of mean distance for all
%               restraints
% rmsd          root mean square error of fitted structure w.r.t. restraints
% network       fitted network
% dxmat         matrix with coordinate displacement vectors for all steps
% diagnostics   .rmsd_vec   vector with restraint rmsd values for each
%                           iteration 
%               .d2         2nd derivative of constraint r.m.s.d.
%               .d2s        smoothed 2nd derivative of constraint r.m.s.d.
%               .converged  iteration number where fit is assumed to have
%                           converged 
%               .drmsd0     initial restraint rmsd
%               .drmsd      final restraint rmsd
%               .timeout    flag, true if fit was stopped by timeout
%               .fit_limit  flag, true if fit reached rmsd limit
%               .maxit      flag, true if fit was stopped by cycle limit
%               .numit      number of iterations
%               .runtime    runtime in seconds
%               
          
% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% G. Jeschke, 2012-2021


% Define defaults
min_fit_limit = ENM_param.tol_model;
maxtime = 3600; % maximum time in seconds

diagnostics.d2 = [];
diagnostics.d2s = [];
diagnostics.converged = 0;
diagnostics.drmsd0 = [];
diagnostics.drmsd = [];
diagnostics.timeout = false;
diagnostics.fit_limit = false;
diagnostics.maxit = 0;
diagnostics.numit = 0;
diagnostics.runtime = 0;

if ~exist('options','var') || isempty(options)
    maxbas = [];
    sat_cyc = 30;
    overfit = false;
else
    if isfield(options,'maxbas')
        maxbas = options.maxbas;
    else
        maxbas = [];
    end
    if isfield(options,'sat_cyc')
        sat_cyc = options.sat_cyc;
    else
        sat_cyc=30;
    end
    if isfield(options,'overfit')
        overfit = options.overfit;
    else
        overfit = false;
    end
end


f = ENM_param.fix_force; % (relative) force constant for restraining local structure
f2 = ENM_param.sec_force; % (relative) force constant for restraining local structure
                        % for second neighbor

[m,~]=size(network0); % m number of residues

dxmat=zeros(m,150);

% set up table of internal restraints for stabilizing local structure
internal=zeros(2*m,4); % initialize table of internal restraints
poi=0;
CaCa_flag = ENM_param.fix_local >= 1;
sec_flag = abs(ENM_param.fix_local) >= 2;
for k = 1:m-2
    cind0 = get_residue_number(CA_model.addresses{k});
    cind1 = get_residue_number(CA_model.addresses{k+1});
    cind2 = get_residue_number(CA_model.addresses{k+2});
    if CaCa_flag && CA_model.chains(k) == CA_model.chains(k+1) && cind1 - cind0 == 1 % consecutive residues
        poi = poi+1;
        internal(poi,1) = k;
        internal(poi,2) = k+1;
        internal(poi,3) = norm(network0(k,:)-network0(k+1,:));
        internal(poi,4) = f;
    end
    if sec_flag && CA_model.chains(k) == CA_model.chains(k+2) && cind2 - cind0 == 2 % residue number difference 2
        poi = poi+1;
        internal(poi,1) = k;
        internal(poi,2) = k + 2;
        internal(poi,3) = norm(network0(k,:) - network0(k+2,:));
        internal(poi,4) = f2;
    end
end
if poi > 0
    internal = internal(1:poi,:);
else
    internal=[];
end

mr = length(ddr); % m number of residues

weights=zeros(mr,1);
for k=1:mr
    weights(k) = 1/ddr(k).sigr;
end
weights = weights.^2;
sc = mr/sum(weights);
weights = sc*weights; % renormalization to conform to Eq. (4) Zheng/Brooks

md=0;

basis = ENM_param.fit_basis;
if basis > 20
    basis=20;
end

fom0 = 0;
weights = ones(size(weights));
for k=1:mr
    vec = ddr(k).xyz2 - ddr(k).xyz1;
    rp = norm(vec);
    fom0 = fom0 + weights(k)*(ddr(k).r-rp)^2/mr;
end

fom0 = sqrt(fom0);
yf = fom0;

rmsd = 0;
for k=1:mr % set up unit vector between restraint residue pairs
    vec = ddr(k).xyz2 - ddr(k).xyz1;
    rp = norm(vec);
    rmsd = rmsd + (ddr(k).r-rp)^2;
end
rmsd = sqrt(rmsd/mr);
rmsd0 = rmsd;
diagnostics.drmsd0 = rmsd0;
fit_limit = 1.9;

numit = 0;
runtime = 0;
dfmax = 0;
ddmax = 0;
rmsd = 1e6;
struct_step = zeros(1,ENM_param.cycles);
change = zeros(1,ENM_param.cycles);
scaling = zeros(1,ENM_param.cycles);
fom_vec = zeros(1,ENM_param.cycles);
rmsd_vec = zeros(1,ENM_param.cycles);
networks = cell(1,ENM_param.cycles);

network00 = network0;
no_stop = 10;
ldfom_vec = zeros(1,ENM_param.cycles);
no_sat = ENM_param.cycles;

% Determine active space extension over cycles
basis_sizes = basis*ones(1,ENM_param.cycles+1);
[ma,~] = size(network0);
basis_offset = ma-6-basis;
bas_incr = linspace(0,basis_offset,ENM_param.cycles+1);
basis_sizes = basis_sizes+round(bas_incr);
if ~isempty(maxbas)
    if maxbas  < basis
        maxbas = basis;
    end
    minbas = basis;
    basis_sizes = maxbas*ones(1,ENM_param.cycles+1);
    basis_sizes(1:sat_cyc) = linspace(minbas,maxbas,sat_cyc);
    basis_sizes = round(basis_sizes);
end

limit_found = false;
ddr0 = ddr;
tic;
while (no_stop || no_sat) && runtime < maxtime && numit < ENM_param.cycles
    if no_stop > 0
        no_stop=no_stop-1;
    else
        no_stop = 0;
    end
    if no_sat > 0
        no_sat = no_sat - 1;
    else
        no_sat = 0;
    end
    [fom,dx,network,ddr_updated,sc] = transition_step(network0,ddr,weights,internal,basis_sizes(numit+1),CA_model,ENM_param);
    bas = 3*numit;
    dxmat(:,bas+1:bas+3) = dx;
    dfom = fom0 - fom;
    numit = numit + 1;
    ldfom_vec(numit) = log(fom0) - log(fom);
    fom_vec(numit) = fom;
    scaling(numit) = sc;
    networks{numit} = network;
    struct_step(numit) = rmsd_superimpose(network,network0);
    change(numit) = rmsd_superimpose(network,network00);
    ddr = ddr_updated;
    if dfom > dfmax 
        dfmax = dfom; 
    end
    if dfom < dfmax/100 && no_sat>10
        no_sat = 10;
    end
    fom0 = fom;
    network0 = network;
    mr = length(ddr);
    rmsd = 0;
    for k = 1:mr % set up unit vector between restraint residue pairs
        vec = ddr(k).xyz2 - ddr(k).xyz1;
        rp = norm(vec);
        rmsd = rmsd + (ddr(k).r-rp)^2;
    end
    rmsd = sqrt(rmsd/mr);
    rmsd_vec(numit) = rmsd;
    if options.interactive
        if numit > 1
            delete(h);
        end
        h = plot(0:numit,[diagnostics.drmsd0 rmsd_vec(1:numit)],'Color',options.color);
        title('Restraint rmsd');
        xlabel('Iteration');
        ylabel(sprintf('restraint rmsd (%s)',char(197)));
    end
    if rmsd < min_fit_limit
        no_stop = 0;
        no_sat = 0;
    end
    ddis = rmsd0 - rmsd;
    if ddis > ddmax, ddmax = ddis; end
    rmsd0 = rmsd;
    axis([0,ENM_param.cycles,0,1.1*max([diagnostics.drmsd0 rmsd_vec])]);
    drawnow;
    if rmsd < fit_limit && ~limit_found
           limit_found = true;
           if ~overfit
              no_stop = 0;
              no_sat = 0;
           end
    end
    runtime = toc;
end

if ~limit_found
    diagnostics.converged = 0;
else
    diagnostics.converged = numit;
end

if numit == ENM_param.cycles
    diagnostics.maxit = true;
else
    diagnostics.maxit = false;
end

best_guess = numit;

if runtime >= maxtime
    diagnostics.timeout = true;
else
    diagnostics.timeout = false;
end
if rmsd <= fit_limit
    diagnostics.fit_limit = true;
else
    diagnostics.fit_limit = false;
end

dxmat=dxmat(:,1:3*numit);

rmsd = 0;
mr = length(ddr);
for k = 1:mr % set up unit vector between restraint residue pairs
    vec = ddr(k).xyz2 - ddr(k).xyz1;
    rp = norm(vec);
    rmsd = rmsd + (ddr(k).r-rp)^2;
end
rmsd = sqrt(rmsd/mr);

if options.interactive
    axis([0,numit,0,1.05*max(yf)]);
    hold on;
    plot([best_guess best_guess],[0,rmsd],'r:');
end
diagnostics.rmsd_vec = rmsd_vec;
improvement = max(rmsd_vec) - min(rmsd_vec);
threshold = min(rmsd_vec) + improvement/3;
der2 = diff(rmsd_vec,2);
diagnostics.d2 = der2;
der2s = der2;
if length(der2) > 1
    der2s(1) = (2*der2(1) + der2(2))/3;
    der2s(end) = (2*der2(end) + der2(end-1))/3;
    for k = 2:length(der2)-1
        der2s(k) = (der2(k-1)+2*der2(k)+der2(k+1))/4;
    end
end
der2s(yf(2:end-1) > threshold) = 0;
diagnostics.d2s = der2s;

thr = fit_limit;
converges = 1;
while converges < length(rmsd_vec) && rmsd_vec(converges) > thr
    converges = converges + 1;
end
if rmsd_vec(converges) == 0
    converges = converges - 1;
end

diagnostics.converged = converges;
plot(converges-1,rmsd_vec(converges),'k^');

diagnostics.drmsd = rmsd_vec(converges);

network=networks{converges-1};

diagnostics.numit = numit;
diagnostics.runtime = runtime;

ddr = update_labels(ddr0,network0,network,CA_model);
rmsd = rmsd_vec(converges-1);

% compute actual error function
fom = 0;
weights = ones(size(weights));
for k=1:mr
    vec = ddr(k).xyz2 - ddr(k).xyz1;
    rp = norm(vec);
    fom = fom+weights(k)*(ddr(k).r-rp)^2/(mr+md);
end
fom = sqrt(fom);


function [fom,dx,network,ddr_updated,sc] = transition_step(network,ddr,weights,internal,basis,CA_model,ENM_param)
% computes one transition step

smallstep = 0.5; % defines an incremental step within the linear regime

mr = length(ddr);

if ENM_param.mass_weighting
    wmat = repmat(sqrt(CA_model.masses'),1,3);
    network = network .* wmat;
end
    
sq_deer=0;
for k = 1:mr
    vec = ddr(k).xyz2 - ddr(k).xyz1;
    rp = norm(vec);
    sq_deer = sq_deer + (ddr(k).r-rp)^2;
end

rmsd = sqrt(sq_deer/mr);
if smallstep < 2*rmsd
    smallstep = rmsd/2;
end

Hessian = get_ANM_Hessian(CA_model,ENM_param);
[u,d] = eig(Hessian);
clear Hessian
lambda = diag(d);
clear d
ut = u(:,7:end);
[~,bas_max] = size(ut);
if basis > bas_max
    basis = bas_max;
end
lambdat = lambda(7:end);
nt = length(lambdat);

network0 = network;
[m,~] = size(network);

if isempty(internal)
    mi=0;
else
    [mi,~]=size(internal);
end

dR = zeros(1,mr+mi);
phases = zeros(nt,1);

for k = 1:mr % set up unit vector between restraint residue pairs
    vec = ddr(k).xyz2 - ddr(k).xyz1;
    rp = norm(vec);
    dR(k) = sqrt(weights(k)/mr)*(ddr(k).r-rp);
end

if mi>0
    for k=1:mi % set up unit vector between local (internal) restraint residue pairs
        vec = network(internal(k,2),:)-network(internal(k,1),:);
        rp = norm(vec);
        dR(k+mr) = sqrt(internal(k,4)/mi)*(internal(k,3)-rp);
    end
end

for j = 1:nt
    du = zeros(1,mr+mi);
    evec = ut(:,j);
    mode = reshape(evec,3,m);
    steps = sqrt(sum(mode.^2,2));
    s = smallstep/max(steps);
    deformed = network0 + s*mode';
    ddr_updated = update_labels(ddr,network0,deformed,CA_model);
    for k=1:mr
        dvec = ddr_updated(k).xyz2 - ddr_updated(k).xyz1;
        drp = norm(dvec);
        vec = ddr(k).xyz2 - ddr(k).xyz1;
        rp = norm(vec);
        du(k) = sqrt(weights(k)/mr)*(drp-rp);
    end
    if mi>0
        for k = 1:mi
            dvec = deformed(internal(k,1),:) - deformed(internal(k,2),:);
            drp = norm(dvec);
            vec = network0(internal(k,2),:) - network0(internal(k,1),:);
            rp = norm(vec);
            du(k+mr) = sqrt(internal(k,4)/mi)*(drp-rp);
        end
    end
    du = du/norm(du);
    dR = dR/norm(dR);
    phases(j) = sum(du.*dR);
end

coeff=sqrt(ones(size(lambdat))./lambdat);

step = ut(:,1:basis)*(phases(1:basis).*coeff(1:basis));
step = reshape(step,3,m);
step = step';
per_residue = sqrt(sum(step.^2,2));
sc = smallstep/max(per_residue);
dx = sc*step;
network = network0 + dx;

ddr_updated = update_labels(ddr,network0,network,CA_model);

% compute actual error function
fom = 0;
weights = ones(size(weights));
for k = 1:mr
    vec = ddr_updated(k).xyz2 - ddr_updated(k).xyz1;
    rp = norm(vec);
    fom = fom + weights(k)*(ddr_updated(k).r-rp)^2/mr;
end
fom = sqrt(fom);

if ENM_param.mass_weighting
    network = network./wmat;
end


function ddr = update_labels(ddr,network0,network,CA_model)
% updates mean spin label cocordinates after a change of C_alpha
% coordinates of the network

[maxnum,~]=size(network);
scarce=0;
% update of label coordinates
for k = 1:length(ddr)
    xyz = ddr(k).xyz1;
    l = ddr(k).index1;
    ind1 = get_residue_number(CA_model.addresses{l});
    local_template = zeros(5,3);
    local_template_0 = zeros(5,3);
    % make a local template to fit rotation and translation
    poi=0;
    for kk = -2:2
        if l+kk > 0 && l+kk <= maxnum % is addressed residue a network point?
            ind2 = get_residue_number(CA_model.addresses{l+kk});
            diff = ind2 - ind1;
            if diff == kk % is addressed residue part of a continuous segment?
                poi = poi+1;
                local_template(poi,:) = network(l+kk,:);
                local_template_0(poi,:) = network0(l+kk,:);
            end
        end
    end
    if poi > 3 % found sufficient number of points to determine local rotation and translation
        [~,~,transmat] = rmsd_superimpose(local_template(1:poi,:),local_template_0(1:poi,:));
        xyz = [xyz 1]; %#ok<AGROW>
        xyz = transmat*xyz';
        xyz = xyz';
        ddr(k).xyz1 = xyz(1:3);
    elseif poi == 3
        [~, ~, transmat] = superimpose_3points(local_template(1:poi,:),local_template_0(1:poi,:));
        xyz = [xyz 1]; %#ok<AGROW>
        xyz = transmat*xyz';
        xyz = xyz';
        ddr(k).xyz1 = xyz(1:3);
    else
        ddr(k).xyz1 = xyz + network(l,:) - network0(l,:); % apply just a translation
        scarce = scarce+1;
    end
    xyz = ddr(k).xyz2;
    l = ddr(k).index2;
    ind1 = get_residue_number(CA_model.addresses{l});
    local_template = zeros(5,3);
    local_template_0 = zeros(5,3);
    % make a local template to fit rotation and translation
    poi=0;
    for kk=-2:2
        if l+kk>0 && l+kk<=maxnum % is addressed residue a network point?
            ind2 = get_residue_number(CA_model.addresses{l+kk});
            diff = ind2 - ind1;
            if diff == kk % is addressed residue part of a continuous segment?
                poi = poi + 1;
                local_template(poi,:) = network(l+kk,:);
                local_template_0(poi,:) = network0(l+kk,:);
            end
        end
    end
    if poi>3 % found sufficient number of points to determine local rotation and translation
        [~,~,transmat] = rmsd_superimpose(local_template(1:poi,:),local_template_0(1:poi,:));
        xyz = [xyz 1]; %#ok<AGROW>
        xyz = transmat*xyz';
        xyz = xyz';
        ddr(k).xyz2 = xyz(1:3);
    elseif poi == 3
        [~, ~, transmat] = superimpose_3points(local_template(1:poi,:),local_template_0(1:poi,:));
        xyz = [xyz 1]; %#ok<AGROW>
        xyz = transmat*xyz';
        xyz=xyz';
        ddr(k).xyz2 = xyz(1:3);
    else
        ddr(k).xyz2 = xyz + network(l,:) - network0(l,:);
        scarce = scarce + 1;
    end
end


function num = get_residue_number(address)

poi = strfind(address,')');
if isempty(poi)
    poi = 0;
end
num = str2double(address(poi+1:end));
