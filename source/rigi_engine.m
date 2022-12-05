function [transforms,diagnostics] = rigi_engine(entity,restraints,options)
% Compute rigid-body arrangements based on distance distribution
% restraints, maximum linker length restraints, cross-links, and clash test 
%
% function [transforms,diagnostics] = RIGI_ENGINE(entity,restraints,options)
%
% Input:
%
% entity        MMMx entity that contains the rigid bodies
% restraints    restraint definition with fields
%               .points     cell(1,rb) of [Nr,3] coordinate vectors of Nr_i 
%                           points in rb rigid bodies that are required
%                           for restraint tests
%               .core       specification of core restraints
%               .lb         lower bounds for core restraints 
%               .ub         upper bounds for core restraints
%               .auxiliary  specification of auxiliary restraints
%               .links      specification of linker length restraints
%               .xlinks     specification of cross-link restraints
%               .r_bodies   definition of rigid bodies in terms of
%                           identifiers of contributing chains
%               .resolution target resolution (Angstroem) for exhaustive
%                           sampling, defaults to 3 Angstroem
%               .p_model    model probability threshold, defaults to 0.5
%               .max_size   maximum size of RBA, defaults to 180 Angstroem
%               .xlff       fraction of cross-links that must be
%                           fulfilled, defaults to 0.3
% 
% options       optional: running options, structure with fields
%               .granule    number of trials in one parallel block,
%                           defaults to 10000
%               .logfid     file identifier for log information, defaults
%                           to Matlab console output
%               .models     maximum number of models, defaults to 20'000
%               .max_time   maximum run time (h), defaults to 48 h
%               .max_trials optional: maximum number of trials, defaults to
%                           no limit, this can also be achieved by setting
%                           it to zero
%               .max_clust  maximum number of RBAs that can be clustered,
%                           default 500000
%
% Output:
%
% transforms    cell(1,M) of cells(1,rb) with affine transformation
%               matrices for rb rigid bodies per model for M models
% diagnostics   struct with diagnostic information on the run with fields
%               .runtime        total run time (s)
%               .trials         number of trials performed
%               .met_err        number of metrization errors encountered
%               .embed_err      number of embedding errors encountered
%               .bound_err      number of bound violations after embedding
%               .aux_fail       auxiliary restraints failures
%               .core_fail      core restraints failures
%               .link_fail      linker length failures
%               .clash_err      rigid-body clashes
%               .xlink_fail     cross-link restraint failures
%               .success        number of successful RBAs
%               .probabilities  model probabilities for successful RBAs
%               .resolution     worst resolution encountered
%               .completed      flag telling whether exhaustive sampling
%                               was completed
%
% G. Jeschke, 2019-2021

% initialize empty output

diagnostics.completed = false;
transforms = [];

if ~exist('options','var') || isempty(options) ...
        || ~isfield(options,'granule') || isempty(options.granule)
    options.granule = 10000;
end

if ~isfield(options,'max_clust') || isempty(options.max_clust)
    options.max_clust = 50000;
end

if ~isfield(restraints,'core')
    restraints.core = [];
end
if ~isfield(restraints,'auxiliary')
    restraints.auxiliary = [];
end
if ~isfield(restraints,'max_trials') || isempty(restraints.max_trials)
    restraints.max_trials = 0;
end

if ~isfield(restraints,'resolution') || isempty(restraints.resolution)
    target_resolution = 3; % resolution with spin labels is not realistically better than 3 Angstroem
else
    target_resolution = restraints.resolution;
end

diagnostics.success = 0;

% settings for the algorithm
maxat = 50000; % maximum number of heavy atoms in rigid body
forgive = 0.8; % scaling factor for van-der-Waals radii in clash tests
min_approach = 5; % minimal approach of two reference points [Angstroem]
max_extension = 180; % maximum distance between any two reference points [Angstroem]
clash_threshold = 1.5*forgive; % a uniform van-der-Waals radius of 1.5 Angstroem is assumed for heavy atoms
clash_fail = 10000; % maximum value of the clash cost function in testing for the unrefined model

if isfield(restraints,'max_size') && ~isempty(restraints.max_size)
    max_extension = restraints.max_size;
end

if isfield(options,'models') && ~isempty(options.models)
    max_models = options.models;
else
    max_models = 20000;
end

if ~isfield(restraints,'p_model') || isempty(restraints.p_model)
    p_model = 0.5;
else
    p_model = restraints.p_model;
end

pthr = exp(-erfinv(p_model)^2);

if isfield(options,'max_time') || ~isempty(options.max_time)
    max_time = options.max_time;
else
    max_time = 48;
end

if isfield(restraints,'xlff') && ~isempty(restraints.xlff)
    xlink_fraction = restraints.xlff;
else
    xlink_fraction = 0.3;
end

if isfield(options,'logfid') && ~isempty(options.logfid)
    logfid = options.logfid;
else
    logfid = 1;
end

% extract full atom coordinates and hevay-atom coordinates for the rigid
% bodies

kc = 0;
all_chains = zeros(1,100);
chain_coor = cell(length(restraints.r_bodies),20);
heavy_coor = cell(1,length(restraints.r_bodies));
for kr = 1:length(restraints.r_bodies)
    coor_r = zeros(maxat,3);
    crpoi = 0;
    for kcc = 1:length(restraints.r_bodies(kr).chains)
        kc = kc + 1;
        all_chains(kc) = restraints.r_bodies(kr).chains(kcc);
        coor = get_coor(entity,sprintf('(%s){1}',all_chains(kc)));
        chain_coor{kr,kcc} = coor;
        coor_r0 = get_coor(entity,sprintf('(%s){1}',all_chains(kc)),true); % heavy-atom coordinates
        [nat,~] = size(coor_r0);
        coor_r(crpoi+1:crpoi+nat,:) = coor_r0;
        crpoi = crpoi + nat;
    end
    heavy_coor{kr} = coor_r(1:crpoi,:);
end

hulls(length(heavy_coor)).vertices = 0;
hulls(length(heavy_coor)).faces = 0;
% convert coordinates to convex hulls
for rba = 1:length(heavy_coor)
    faces = convhulln(heavy_coor{rba});
    [na,~] = size(faces);
    [kpa,ar] = reducepatch(faces,heavy_coor{rba},na);
    hulls(rba).vertices = double(ar);
    hulls(rba).faces = kpa;
end

% number of rigid bodies
rbnum = length(restraints.r_bodies);

% augment lower and upper bounds
for k1 = 1:3*rbnum-1
    for k2 = k1+1:3*rbnum
        if restraints.lb(k1,k2) < min_approach
            restraints.lb(k1,k2) = min_approach;
            restraints.lb(k2,k1) = min_approach;
        end
        if restraints.ub(k1,k2) < 0.1 % the unset bounds are zero
            restraints.ub(k1,k2) = max_extension;
            restraints.ub(k2,k1) = max_extension;
        end
    end
end

[trials,res,trial_pattern] = get_restraint_resolution(restraints.lb,restraints.ub,target_resolution);
if options.max_trials > 0
    while trials > options.max_trials
        target_resolution = target_resolution + 0.1;
        [trials,res,trial_pattern] = get_restraint_resolution(restraints.lb,restraints.ub,target_resolution);
    end
end
intervals = trial_pattern(:,3);
digitbase = trial_pattern(:,4);
irbr = trial_pattern(:,1:2);

[lb,ub,err]=triangle(restraints.lb,restraints.ub);

switch err
    case 0
        fprintf(logfid,'Successful bound smoothing with experimental restraints.\n');
    case 1
        fprintf(logfid,'ERROR: Some distance restraints are inconsistent.\n');
        diagnostics.success = -1;
        return
    otherwise
        fprintf(logfid,'Unspecified error in bound smoothing.\n');
        diagnostics.success = -1;
        return
end

% initalize failure counters
met_err = 0;
embed_err = 0;
bound_err = 0;
clash_err = 0;
aux_fail = 0;
link_fail = 0;
xlink_fail = 0;
core_fail = 0;
success = 0;


% the following reduces communication overhead for the parfor loop by not
% supplying the whole restraints variable
rb = restraints.r_bodies;
points = restraints.points;
core = restraints.core;
auxiliary = restraints.auxiliary;
links = restraints.links;
xlinks = restraints.xlinks;
[naux,~] = size(auxiliary);
[ncore,~] = size(core);

transmats = cell(1,3);

% initialize counters and diagnostic variables
runtime = 0;
worst_res = 0;
parblocks = 0; 
bask = parblocks*options.granule;

probabilities = zeros(1,options.max_clust);
transforms = zeros(options.max_clust,6*length(rb));
all_ref_points = cell(1,options.max_clust);
ref_points = zeros(3*length(rb),3);
for kb = 1:length(rb)
    bpoints = points{kb};
    bas = 3*(kb-1);
    ref_points(bas+1:bas+3,:) = bpoints(1:3,:);
end

tic,

while bask < trials && runtime <= 3600*max_time && success < options.max_clust
    % initialize diagnostic variable vectors for one parallel block
    merr_vec = zeros(1,options.granule);
    eerr_vec = zeros(1,options.granule);
    berr_vec = zeros(1,options.granule);
    aerr_vec = zeros(1,options.granule);
    rerr_vec = zeros(1,options.granule);
    lerr_vec = zeros(1,options.granule);
    cerr_vec = zeros(1,options.granule);
    xerr_vec = zeros(1,options.granule);
    model_prob = zeros(1,options.granule);
   
    tmats = cell(options.granule,1);
    parblocks = parblocks + 1;
    res_vec = res*ones(1,options.granule);
    % profile on
    parfor kt = 1:options.granule % ### edit to 'for' if you want to debug
        fulfill = true;
        if bask + kt <= prod(intervals)
            fractions = get_restraint_fractions(bask+kt,intervals,digitbase);
            [dmatr,err,cres] = metrize_exhaustive(lb,ub,irbr,fractions,intervals);
        else
            err = 1; % report a metrization error for surplus trials
        end
        if err == 1% metrization failed (restraints inconsistent), next trial, increment error counter
            merr_vec(kt) = 1;
        else
            if cres > res_vec(kt)
                res_vec(kt) = cres;
            end
            coor0 = dmat2coor(dmatr); % embed distance matrix to Cartesian space
            if isempty(coor0)
                eerr_vec(kt) = 1;
            else
                [coor1,err] = bound_refiner(coor0,lb,ub);
                if err > 0
                    berr_vec(kt) = 1;
                else
                    t_labels = cell(1,length(rb));
                    t_points = cell(1,length(rb));
                    atransmat = zeros(4*length(rb),4);
                    for kr = 1:length(rb)
                        baspoi = 3*(kr-1);
                        baspoi4 = 4*(kr-1);
                        [~,lab_t,transmat] = superimpose_3points(coor1(baspoi+1:baspoi+3,:),rb(kr).ref(1:3,2:4));
                        t_labels{kr} = lab_t;
                        t_points{kr} = affine_coor_set(points{kr},transmat); %#ok<PFBNS>
                        atransmat(baspoi4+1:baspoi4+4,:) = transmat;
                    end
                    tmats{kt} = atransmat;
                    paux = 1;
                    for kaux = 1:naux
                        acoor1 = t_points{auxiliary(kaux,1)}; %#ok<PFBNS>
                        acoor2 = t_points{auxiliary(kaux,3)};
                        coor1b = acoor1(auxiliary(kaux,2),:);
                        coor2b = acoor2(auxiliary(kaux,4),:);
                        rsim = norm(coor1b-coor2b);
                        prob = prob_Gaussian(rsim,auxiliary(kaux,5),sqrt(auxiliary(kaux,6)^2+cres^2)); % include resolution into uncertainty
                        paux = paux*prob;
                    end
                    if paux < pthr^naux
                        aerr_vec(kt) = 1;
                        fulfill = false;
                    end
                    % check how well label-to-label distances are reproduced
                    plabels = 1;
                    if fulfill
                        for kcore = 1:ncore
                            acoor1 = t_points{core(kcore,1)}; %#ok<PFBNS>
                            acoor2 = t_points{core(kcore,3)};
                            coor1b = acoor1(core(kcore,2),:);
                            coor2b = acoor2(core(kcore,4),:);
                            rsim = norm(coor1b-coor2b);
                            prob = prob_Gaussian(rsim,core(kcore,5),sqrt(core(kcore,6)^2+cres^2)); % include resolution into uncertainty
                            plabels = plabels*prob;
                        end
                        if plabels < pthr^ncore
                            fulfill = false;
                            rerr_vec(kt) = 1;
                        end
                    end
                    % check for linker lengths
                    if fulfill
                        if ~isempty(links(1).maxr) %#ok<PFBNS>
                            for kl = 1:length(links)
                                if ~isempty(links(kl).ref_indices)
                                    acoor1 = t_points{links(kl).ref_indices(1)};
                                    acoor2 = t_points{links(kl).ref_indices(3)};
                                    coor1b = acoor1(links(kl).ref_indices(2),:);
                                    coor2b = acoor2(links(kl).ref_indices(4),:);
                                    rlink = norm(coor2b - coor1b);
                                    if rlink > links(kl).maxr + cres % include resolution into uncertainty
                                        fulfill = false;
                                    end
                                end
                            end
                            if ~fulfill
                                lerr_vec(kt) = 1;
                            end
                        end
                    end
                    % check for rigid-body clashes
                    if fulfill
                        clashed = false;
                        for kr1 = 1:length(rb)-1
                            baspoi4 = 4*(kr1-1);
                            atransmat = tmats{kt};
                            transmat = atransmat(baspoi4+1:baspoi4+4,:);
                            hc1 = single(affine_coor_set(heavy_coor{kr1},transmat)); %#ok<PFBNS>
                            hc1h = affine_coor_set(hulls(kr1).vertices,transmat); %#ok<PFBNS>
                            for kr2 = kr1+1:length(rb)
                                baspoi4b = 4*(kr2-1);
                                transmatb = atransmat(baspoi4b+1:baspoi4b+4,:);
                                hc2 = single(affine_coor_set(heavy_coor{kr2},transmatb));
                                hc2h = affine_coor_set(hulls(kr2).vertices,transmatb);
                                cost = clash_cost_super_fast(hc1h,hulls(kr1).faces,hc2h,hulls(kr2).faces,hc1,hc2,clash_threshold,-1);
                                if cost > clash_fail
                                    clashed = true;
                                end
                            end
                        end
                        if clashed
                            cerr_vec(kt) = 1;
                            fulfill = false;
                        end
                    end
                    % check for cross-link fulfillment
                    if fulfill
                        xld = zeros(length(xlinks),1);
                        xlf = 0;
                        for kl = 1:length(xlinks)
                            acoor1 = t_points{xlinks(kl).ref_indices(1)};
                            acoor2 = t_points{xlinks(kl).ref_indices(3)};
                            coor1b = acoor1(xlinks(kl).ref_indices(2),:);
                            coor2b = acoor2(xlinks(kl).ref_indices(4),:);
                            rlink = norm(coor2b - coor1b);
                            xld(kl) = rlink;
                            if rlink > xlinks(kl).maxr + cres % include resolution into uncertainty
                                xlf = xlf + 1;
                            end
                        end
                        if (length(xlinks)-xlf)/length(xlinks) < xlink_fraction
                            xerr_vec(kt) = 1;
                            fulfill = false;
                        end
                    end
                    if fulfill
                        atransmat = tmats{kt};
                        % the first refinement is performed without
                        % testing for clashes, as this is much faster
                        atransmat0 = atransmat;
                        [atransmat,errvec,mprob] = ...
                            refine_rba(rb,atransmat0,points,pthr,naux,auxiliary,ncore,core,links,1e6,heavy_coor,100*xlink_fraction,xlinks);
                        if ~sum(errvec)
                            [atransmat,errvec,mprob] = ...
                                refine_rba_fast(rb,atransmat,points,pthr,naux,auxiliary,ncore,core,links,clash_threshold,heavy_coor,100*xlink_fraction,xlinks);
                        end
                        tmats{kt} = atransmat;
                        model_prob(kt) = mprob;
                        aerr_vec(kt) = errvec(1);
                        rerr_vec(kt) = errvec(2);
                        lerr_vec(kt) = errvec(3);
                        cerr_vec(kt) = errvec(4);
                        xerr_vec(kt) = errvec(5);
                    end
                end
            end
        end
    end
    fprintf(logfid,'Parallel block %i completed.\n',parblocks);
    % update the failure counters
    met_err = met_err + sum(merr_vec);
    embed_err = embed_err + sum(eerr_vec);
    bound_err = bound_err + sum(berr_vec);
    aux_fail = aux_fail + sum(aerr_vec);
    core_fail = core_fail + sum(rerr_vec);
    link_fail = link_fail + sum(lerr_vec);
    clash_err = clash_err + sum(cerr_vec);
    xlink_fail = xlink_fail + sum(xerr_vec);
    % generate the vector of failure flags for all parallel trials
    fail_vec = merr_vec+eerr_vec+berr_vec+aerr_vec+rerr_vec+lerr_vec+cerr_vec+xerr_vec;
    % update the worts reolution encountered
    if worst_res < max(res_vec)
        worst_res = max(res_vec);
    end
    % process all results from this parallel block
    for k = bask+1:bask+options.granule
        if fail_vec(k-bask) == 0 % this is a successfully refined RBA
            atransmat = tmats{k-bask};
            tref_points = zeros(3*length(rb),3);
            for kr = 1:length(rb)
                baspoi4 = 4*(kr-1);
                transmats{kr} = atransmat(baspoi4+1:baspoi4+4,:);
                baspoi3 = 3*(kr-1);
                bpoints = ref_points(1+baspoi3:3+baspoi3,:);
                tref_points(1+baspoi3:3+baspoi3,:) = affine_coor_set(bpoints,transmats{kr});
            end
            success = success + 1;
            if success <= options.max_clust
                probabilities(success) = model_prob(k-bask)^(1/(naux+ncore));
                transvecs = zeros(1,6*length(rb));
                all_ref_points{success} = tref_points;
                for kr = 1:length(rb)
                    bas = 6*(kr-1);
                    [trans,euler] = affine2transrot(transmats{kr});
                    transvecs(bas+1:bas+3) = trans;
                    transvecs(bas+4:bas+6) = euler;
                end
                transforms(success,:) = transvecs;
            end
        end
    end
    bask = bask + options.granule;
    runtime=toc;
    dmg_err = met_err + embed_err + bound_err;
    %  update runtime information
    left_time1 = 3600*max_time - runtime;
    left_time2 = runtime*(trials-bask)/bask;
    left_time3 = runtime*(options.max_clust-success)/success;
    if isnan(left_time3)
        left_time3 = max_time;
    end
    % Report current state to logfile
    left_time = min([left_time1 left_time2 left_time3]);
    done_trials = bask;
    done_trials = done_trials - dmg_err;
    hours = floor(left_time/3600);
    minutes = round((left_time-3600*floor(left_time/3600))/60);
    fprintf(logfid,'%i h %i min estimated run time to completion.\n',hours,minutes);
    fprintf(logfid,'%6.2f%% completed.\n',100*bask/trials);
    fprintf(logfid,'%i feasible rigid-body arrangements found.\n',success);
    fprintf(logfid,'%12.1f s computation time per model.\n',runtime/success);
    fprintf(logfid,'%6.2f%% distance geometry failures.\n',100*dmg_err/bask);
    fprintf(logfid,'%6.2f%% auxiliary restraint failuers.\n',100*aux_fail/done_trials);
    fprintf(logfid,'%6.2f%% core restraint failures.\n',100*core_fail/done_trials);
    fprintf(logfid,'%6.2f%% linker length failures.\n',100*link_fail/done_trials);
    fprintf(logfid,'%6.2f%% rigid-body clashes.\n',100*clash_err/done_trials);
    fprintf(logfid,'%6.2f%% cross-link failures.\n',100*xlink_fail/done_trials);
end
toc,

solutions = success;

transforms = transforms(1:success,:);
probabilities = probabilities(1:success);

cluster_resolution = NaN;
cluster_time = 0;
if success > max_models
    tic,
    cluster_resolution = 0;
    % hierarchical clustering
    % make distance matrix between 
    dmat = zeros(success);
    for k1 = 1:success-1
        coor1 = all_ref_points{k1};
        for k2 = k1+1:success
            coor2 = all_ref_points{k2};
            rms = rmsd_superimpose(coor1,coor2);
            dmat(k1,k2) = rms;
            dmat(k2,k1) = rms;
        end
    end
    % compute linkage
    Z = linkage(dmat,'average');
    % cluster
    assignment = cluster(Z,'maxclust',max_models);
    % select representative model for each cluster
    selected_models = zeros(max_models,1);
    all_indices = 1:success;
    for km = 1:max_models
        % determine indices corresponding to cluster km
        indices = all_indices(assignment==km);
        % extract distance matrix for this cluster
        dmat_cluster = dmat(indices,indices);
        % compute mean square sum for each model in cluster with respect to
        % all other models
        msq_models = sum(dmat_cluster.^2);
        % find minimum mean square deviation and index of this model in
        % cluster
        [~,index] = min(msq_models);
        selected_models(km) = indices(index);
        % maximum rooot mean square distance of another model in this
        % cluster from selected model
        max_rmsd = max(dmat_cluster(index,:));
        if max_rmsd > cluster_resolution
            cluster_resolution = max_rmsd;
        end
    end
    transforms = transforms(selected_models,:);
    probabilities = probabilities(selected_models);
    success = max_models;
    cluster_time = toc;
end

diagnostics.runtime = runtime;
diagnostics.trials = bask;
diagnostics.met_err = met_err;
diagnostics.embed_err = embed_err;
diagnostics.bound_err = bound_err;
diagnostics.aux_fail = aux_fail;
diagnostics.core_fail = core_fail;
diagnostics.link_fail = link_fail;
diagnostics.clash_err = clash_err;
diagnostics.xlink_fail = xlink_fail;
diagnostics.success = success;
diagnostics.solutions = solutions;
diagnostics.probabilities = probabilities;
diagnostics.resolution = worst_res;
diagnostics.cluster_resolution = cluster_resolution;
diagnostics.cluster_time = cluster_time;
diagnostics.completed = bask >= trials;

