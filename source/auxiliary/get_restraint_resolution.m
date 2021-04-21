function [ntrials,res,trial_pattern,minlb,maxub] = get_restraint_resolution(lb,ub,restarget)

if ~exist('restarget','var') || isempty(restarget)
    restarget = 1e6;
end

ntrials = 1;
res = 0;
minlb = 1e6;
maxub = 0;

[n,~] = size(lb);
rb = n/3; % number of rigid bodies, assuming three reference points per rigid body
irbr = 9*rb*(rb-1)/2; % number of inter-rigid-body restraints

trial_pattern = zeros(irbr,4);

digitbase = 1;
for k = 1:irbr
    [lb,ub,err] = triangle_bound_smoothing(lb,ub,1);
    if err
        res = [];
        ntrials = [];
        return
    end
    range = ub - lb;
    range(range == 0) = 10000;
    [mi,poi1] = min(range);
    [mi,poi2] = min(mi);
    poi1 = poi1(poi2);
    ctrials = 1;
    while mi/(ctrials+1) > restarget
        ctrials = ctrials + 1;
    end
    if mi/(ctrials+1) > res
        res = mi/(ctrials+1);
    end
    ntrials = ntrials * ctrials;
    if ub(poi1,poi2) > maxub
        maxub = ub(poi1,poi2);
    end
    if lb(poi1,poi2) < minlb
        minlb = lb(poi1,poi2);
    end
    lb(poi1,poi2) = lb(poi1,poi2) + range(poi1,poi2)/2;
    lb(poi2,poi1) = lb(poi1,poi2);
    ub(poi1,poi2) = lb(poi1,poi2);
    ub(poi2,poi1) = lb(poi1,poi2);
    trial_pattern(k,:) = [poi1,poi2,ctrials,digitbase];
    digitbase = digitbase * ctrials;
end