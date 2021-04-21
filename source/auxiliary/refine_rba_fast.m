function [atransmat,errvec,model_prob,xld,cost,costs,aux_fulfill] = refine_rba_fast(rb,atransmat,points,...
    pthr,naux,auxiliary,ncore,core,links,clash_threshold,heavy_coor,xlink_percentage,xlinks,...
    stemlibs,combination)

global v_first_fit
global rba_refined

v_first_fit = [];
rba_refined = false;

clash_threshold = single(clash_threshold);

verbose = false;
detailed_verbose = false;
full_verbose = false;
preliminary_verbose = false;
do_echo = false;

clash_fail = 10; % value of the clash cost function were clashes are considered as too strong 

if ~exist('stemlibs','var')
    stemlibs = {};
end

if ~exist('combination','var')
    combination = [];
end

if ~isempty(combination)
    % augment heavy_coor with stemlib coordinates if required
    for klib = 1:length(stemlibs)
        rba = stemlibs{klib}.rba;
        sl_coor = stemlibs{klib}.chains{combination(klib)}.xyz{1};
        elements = stemlibs{klib}.chains{combination(klib)}.isotopes;
        sl_coor = sl_coor(elements(:,1) ~= 1,:);
%         baspoi4 = 4*(rba-1);
%         transmat = atransmat(baspoi4+1:baspoi4+4,:);
%         sl_coor = single(affine_coor_set(sl_coor,transmat));
        coor0 = heavy_coor{rba};
        heavy_coor{rba} = [coor0; sl_coor];
    end
end

% convert coordinates to convex hulls
for rba = 1:length(heavy_coor)
    faces = convhulln(double(heavy_coor{rba}));
    [na,~] = size(faces);
    [kpa,ar] = reducepatch(faces,double(heavy_coor{rba}),na);
    hulls(rba).vertices = double(ar);
    hulls(rba).faces = kpa;
end

sl_atoms = 0;
for klib = 1:length(stemlibs)
    sl_coor = stemlibs{klib}.chains{1}.xyz{1};
    [msla,~] = size(sl_coor);
    sl_atoms = sl_atoms + msla;
end

errvec = zeros(1,5);
t_labels = cell(1,length(rb));
t_points = cell(1,length(rb));

xld = [];

v0 = zeros(1,6*(length(rb)-1)); % vector of Euler angles and translation vectors
% fprintf(1,'%i rigid bodies\n',length(rb));
% fprintf(1,'Size of the transformation matrix is (%i,%i)\n',size(atransmat));
for kr = 2:length(rb)
    baspoi4 = 4*(kr-1);
    baspoi6 = 6*(kr-2);
    transmat = atransmat(baspoi4+1:baspoi4+4,:);
    [trans,euler] = affine2transrot(transmat);
    v0(baspoi6+1:baspoi6+3) = trans/10; % go to nm, better balance of x values
    v0(baspoi6+4:baspoi6+6) = euler;
end

options = optimset('Display','none','TolFun',1,'TolX',0.03,'MaxFunEvals',10000);

% v1 = fminsearch(@rba_cost_fct,v0,options,atransmat,rb,points,naux,auxiliary,ncore,core,links,clash_threshold,hulls);

% lb = v0 - 0.6*ones(size(v0));
% ub = v0 + 0.6*ones(size(v0));

if ~isempty(combination) && preliminary_verbose
    [~,costs] = rba_cost_fct(v0,atransmat,rb,points,naux,auxiliary,ncore,core,links,clash_threshold,heavy_coor,hulls,pthr,clash_fail,false,v0);
    fprintf(1,'Optimization (%i,%i,%i) started with aux. %5.2f; core %5.2f; links %5.2f; clash %5.2f; stem %5.2f\n',combination,costs.aux,costs.core,costs.link,costs.clash,costs.stem);
end

tstart = tic;
if ~isempty(combination)
    [v1,cost] = fminsearch(@rba_cost_fct,v0,options,atransmat,rb,points,naux,auxiliary,ncore,core,links,clash_threshold,heavy_coor,hulls,pthr,clash_fail,do_echo,v0);
else
    [v1,cost] = fminsearch(@rba_cost_fct,v0,options,atransmat,rb,points,naux,auxiliary,ncore,core,links,clash_threshold,heavy_coor,hulls,pthr,clash_fail,do_echo,v0);
end
telapsed = toc(tstart);

% rba_refined = false;

if rba_refined
    v1 = v_first_fit;
end

Dv = v1 - v0;
if length(Dv) == 12
    trans = [Dv(1:3) Dv(7:9)];
    maxtrans = max(abs(trans));
    rot1 = [Dv(4:6) Dv(10:12)];
    rot2 = 2*pi-rot1;
    maxrot = 0;
    for k = 1:6
        ma1 = min(abs([rot1(k) rot2(k)]));
        if ma1 > maxrot
            maxrot = ma1;
        end
    end
else
    maxtrans = NaN;
    maxrot = NaN;
end
[cost,costs] = rba_cost_fct(v1,atransmat,rb,points,naux,auxiliary,ncore,core,links,clash_threshold,heavy_coor,hulls,pthr,clash_fail,false,v0);
if ~isempty(combination)
    if (clash_threshold < 2 && full_verbose) || preliminary_verbose
        fprintf(1,'Optimization (%i,%i,%i) converged to cost %6.1f with max. trans. %5.2f nm and max. rot. %5.2f degree (%5.1f s)\n',combination,cost,maxtrans,180*maxrot/pi,telapsed);
        fprintf(1,'aux. %5.2f; core %5.2f; links %5.2f; clash %5.2f; stem %5.2f\n',costs.aux,costs.core,costs.link,costs.clash,costs.stem);
    end
end

% if clash_threshold < 2
%     disp('Aber hallo!');
% end


model_prob = 0;

for kr = 1:length(rb)
    baspoi4 = 4*(kr-1);
    baspoi6 = 6*(kr-2);
    if kr > 1
        trans = 10*v1(baspoi6+1:baspoi6+3); % reconvert to Angstroem
        euler = v1(baspoi6+4:baspoi6+6);
        transmat = transrot2affine(trans,euler);
    else
        transmat = atransmat(baspoi4+1:baspoi4+4,:);
    end
    t_labels{kr} = affine_coor_set(rb(kr).ref(1:3,2:4),transmat);
    t_points{kr} = affine_coor_set(points{kr},transmat);
    atransmat(baspoi4+1:baspoi4+4,:) = transmat;
end

paux = 1;
aux_fulfill = zeros(1,naux);
for kaux = 1:naux
    acoor1 = t_points{auxiliary(kaux,1)};
    acoor2 = t_points{auxiliary(kaux,3)};
    coor1b = acoor1(auxiliary(kaux,2),:);
    coor2b = acoor2(auxiliary(kaux,4),:);
    rsim = norm(coor1b-coor2b);
    prob = prob_Gaussian(rsim,auxiliary(kaux,5),auxiliary(kaux,6)); 
    if detailed_verbose 
        if clash_threshold >=2
            fprintf(1,'PRE: ');
        else
            fprintf(1,'FIN: ');
        end
    end
    aux_fulfill(kaux) = prob;
    if detailed_verbose
        fprintf(1,'Auxiliary restraint %i at %4.1f ? for %4.1f ? (%4.1f ?); Probability %4.2f\n',kaux,rsim,auxiliary(kaux,5),auxiliary(kaux,6),prob);
    end
    paux = paux*prob;
end
% fprintf(1,'Auxiliary restraint fulfilment: %8.5f\n',paux/pthr^naux);
if paux < pthr^naux && ~rba_refined
    if (clash_threshold < 2 && verbose) || preliminary_verbose
        fprintf(2,'Auxiliary restraint fulfilment: %8.5f\n',paux/pthr^naux)
    end
    errvec(1) = 1;
    return
elseif (clash_threshold < 2 && verbose) || preliminary_verbose
    fprintf(1,'Auxiliary restraint fulfilment: %8.5f\n',paux/pthr^naux)
end

% check how well label-to-label distances are reproduced
plabels = 1;
for kcore = 1:ncore
    acoor1 = t_points{core(kcore,1)};
    acoor2 = t_points{core(kcore,3)};
    coor1b = acoor1(core(kcore,2),:);
    coor2b = acoor2(core(kcore,4),:);
    rsim = norm(coor1b-coor2b);
    prob = prob_Gaussian(rsim,core(kcore,5),core(kcore,6));
    plabels = plabels*prob;
end
% fprintf(1,'Core restraint fulfilment: %8.5f\n',plabels/pthr^ncore);
if plabels < pthr^ncore && ~rba_refined
    if (clash_threshold < 2 && verbose) || preliminary_verbose
        fprintf(2,'Core restraint fulfilment: %8.5f\n',plabels/pthr^ncore);
    end
    errvec(2) = 1;
    return;
end

model_prob = paux^(1/naux)*plabels^(1/ncore);
% check for linker lengths
if ~isempty(links(1).maxr)
    fulfill = true;
    for kl = 1:length(links)
        if ~isempty(links(kl).ref_indices)
            acoor1 = t_points{links(kl).ref_indices(1)};
            acoor2 = t_points{links(kl).ref_indices(3)};
            coor1b = acoor1(links(kl).ref_indices(2),:);
            coor2b = acoor2(links(kl).ref_indices(4),:);
            rlink = norm(coor2b - coor1b);
            if rlink > links(kl).maxr 
                fulfill = false;
            end
        end
    end
    if ~fulfill && ~rba_refined
        if (clash_threshold < 2 && verbose) || preliminary_verbose
            fprintf(2,'Linker lengths not fulfilled.\n');
        end
        errvec(3) = 1;
        return
    end
end

xld = zeros(length(xlinks),1);
xlf = 0;
for kl = 1:length(xlinks)
    acoor1 = t_points{xlinks(kl).ref_indices(1)};
    acoor2 = t_points{xlinks(kl).ref_indices(3)};
    coor1b = acoor1(xlinks(kl).ref_indices(2),:);
    coor2b = acoor2(xlinks(kl).ref_indices(4),:);
    rlink = norm(coor2b - coor1b);
    xld(kl) = rlink;
    if rlink > xlinks(kl).maxr
        xlf = xlf + 1;
    end
end
% fprintf(1,'Crosslink restraint fulfilment: %6.3f\n',(100*(length(xlinks)-xlf)/length(xlinks))/xlink_percentage);
if 100*(length(xlinks)-xlf)/length(xlinks) < xlink_percentage
    errvec(5) = 1;
    return
end

% check for rigid-body clashes
clashed = false;
max_cost = 0;
for kr1 = 1:length(rb)-1
    baspoi4 = 4*(kr1-1);
    transmat = atransmat(baspoi4+1:baspoi4+4,:);
    hc1 = single(affine_coor_set(heavy_coor{kr1},transmat));
    hc1h = single(affine_coor_set(hulls(kr1).vertices,transmat));
    for kr2 = kr1+1:length(rb)
        baspoi4b = 4*(kr2-1);
        transmatb = atransmat(baspoi4b+1:baspoi4b+4,:);
        hc2 = single(affine_coor_set(heavy_coor{kr2},transmatb));
        hc2h = single(affine_coor_set(hulls(kr2).vertices,transmatb));
        rcost = clash_cost_super_fast(hc1h,hulls(kr1).faces,hc2h,hulls(kr2).faces,hc1,hc2,clash_threshold);
        if rcost > max_cost
            max_cost = rcost;
        end
        if rcost > clash_fail
            clashed = true;
        end
    end
end
if clashed && ~rba_refined
    if clash_threshold < 2 && verbose
        fprintf(2,'Core rigid body clashes persist at cost %8.1f (threshold %3.1f) with limit %4.0f.\n',max_cost,clash_threshold,clash_fail);
    end
    errvec(4) = 1;
elseif clash_threshold < 1e4 && verbose
    fprintf(1,'Core clashes reduced to cost %8.1f (threshold %3.1f) with limit %4.0f.\n',max_cost,clash_threshold,clash_fail);
end


function [cost,costs] = rba_cost_fct(v,atransmat,rb,points,naux,auxiliary,ncore,core,links,clash_threshold,heavy_coor,hulls,pthr,clash_fail,echo,v0)

global v_first_fit
global rba_refined

link_weight = 1;

fit_success = true;

if ~exist('echo','var')
    echo = false;
end
if echo
    fprintf(1,'Trafo 1: (%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f)\n ',v(1:6)-v0(1:6));
    fprintf(1,'Trafo 2: (%5.2f,%5.2f,%5.2f,%5.2f,%5.2f,%5.2f)\n ',v(7:12)-v0(7:12));
end
cost = 0;

t_labels = cell(1,length(rb));
t_points = cell(1,length(rb));

for kr = 1:length(rb)
    baspoi4 = 4*(kr-1);
    baspoi6 = 6*(kr-2);
    if kr > 1
        trans = 10*v(baspoi6+1:baspoi6+3); % reconvert to Angstroem
        euler = v(baspoi6+4:baspoi6+6);
        transmat = transrot2affine(trans,euler);
    else
        transmat = atransmat(baspoi4+1:baspoi4+4,:);
    end
    t_labels{kr} = affine_coor_set(rb(kr).ref(1:3,2:4),transmat);
    t_points{kr} = affine_coor_set(points{kr},transmat);
    atransmat(baspoi4+1:baspoi4+4,:) = transmat;
end

% auxiliary restraint chi^2
cost_aux = 0;
for kaux = 1:naux
    acoor1 = t_points{auxiliary(kaux,1)};
    acoor2 = t_points{auxiliary(kaux,3)};
    coor1b = acoor1(auxiliary(kaux,2),:);
    coor2b = acoor2(auxiliary(kaux,4),:);
    rsim = norm(coor1b-coor2b);
    cost_aux = cost_aux + ((rsim-auxiliary(kaux,5))/auxiliary(kaux,6))^2; 
end
if naux > 0
    costs.aux = cost_aux/(-naux*log(pthr));
else
    costs.aux = 0;
end
cost = cost + costs.aux;
if costs.aux > 1
    fit_success = false;
end
if echo
    fprintf(1,'Auxiliary : %5.2f\n',costs.aux);
end

% core restraint chi^2
cost_core = 0;
for kcore = 1:ncore
    acoor1 = t_points{core(kcore,1)};
    acoor2 = t_points{core(kcore,3)};
    coor1b = acoor1(core(kcore,2),:);
    coor2b = acoor2(core(kcore,4),:);
    rsim = norm(coor1b-coor2b);
    cost_core = cost_core + ((rsim-core(kcore,5))/core(kcore,6))^2; 
end
if ncore > 0
    costs.core = cost_core/(-ncore*log(pthr));
else
    costs.core = 0;
end
cost = cost + costs.core;
if costs.core > 1
    fit_success = false;
end
if echo
    fprintf(1,'Core      : %5.2f\n',costs.core);
end

% check for linker lengths
cost_link = 0;
if ~isempty(links(1).maxr)
    for kl = 1:length(links)
        if ~isempty(links(kl).ref_indices)
            acoor1 = t_points{links(kl).ref_indices(1)};
            acoor2 = t_points{links(kl).ref_indices(3)};
            coor1b = acoor1(links(kl).ref_indices(2),:);
            coor2b = acoor2(links(kl).ref_indices(4),:);
            rlink = norm(coor2b - coor1b);
            % fprintf(1,'Link(%i), Combi(%i,%i,%i): %4.2f ?\n',kl,combination,rlink);
            if rlink > 0.99*links(kl).maxr 
                cost_link = cost_link + link_weight*(rlink-0.99*links(kl).maxr)^2;
            end
        end
    end
end
costs.link = cost_link;
cost = cost + cost_link;
if cost_link > 1
    fit_success = false;
end
if echo
    fprintf(1,'Links     : %5.2f\n',cost_link);
end

costs.clash = 0;
costs.stem = 0;
if clash_threshold < 1e4
    hcs = cell(1,length(rb));
    hcsh = cell(1,length(rb));
    mcr = zeros(1,length(rb));
    cost_clash = 0;
    for kr = 1:length(rb)
        baspoi4b = 4*(kr-1);
        transmatb = atransmat(baspoi4b+1:baspoi4b+4,:);
        hcs{kr} = single(affine_coor_set(heavy_coor{kr},transmatb));
        hcsh{kr} = single(affine_coor_set(hulls(kr).vertices,transmatb));
        [m,~] = size(hcs{kr});
        mcr(kr) = m;
    end
    for kr1 = 1:length(rb)-1
        hc1h = hcsh{kr1};
        hc1 = hcs{kr1};
        for kr2 = kr1+1:length(rb)
            hc2h = hcsh{kr2};
            hc2 = hcs{kr2};
            rcost = clash_cost_super_fast(hc1h,hulls(kr1).faces,hc2h,hulls(kr2).faces,hc1,hc2,clash_threshold);
            if rcost > cost_clash
                cost_clash = rcost;
            end
        end
    end
    costs.clash = cost_clash/clash_fail;
    cost = cost + costs.clash;
    if costs.clash > 1
        fit_success = false;
    end
    if echo
        fprintf(1,'Core clash: %5.2f\n',costs.clash);
    end
end
if echo
    fprintf(2,'Total cost: %5.2f\n',cost);
end

if fit_success && ~rba_refined
    rba_refined = true;
    v_first_fit = v;
end