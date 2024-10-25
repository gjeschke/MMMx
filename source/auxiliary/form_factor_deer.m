function ff = form_factor_deer(r0,Pr0,t0,Pake_base)

Pr0 = Pr0/sum(Pr0);
dt = Pake_base.t(2)-Pake_base.t(1);
dt0 = t0(2)-t0(1);
stretch = (dt0/dt)^(1/3);
try
    distr = get_std_distr(r0/stretch,Pr0,Pake_base.r);
catch
    ff = [];
    return
end
distr=0.01*distr/sum(distr);
rsig = distr*Pake_base.kernel; % multiply G(r) vector with kernel matrix
ff0 = exp(rsig); % eqn [10]
if sum(isnan(ff0))
    ff = ones(size(t0));
else
    ff = interp1(Pake_base.t*stretch^3,ff0,t0,'pchip',0.99);
    ff = ff - 0.99;
    ff = ff/ff(1);
end

function distr = get_std_distr(r0,p0,r)
%
% Computes distance distribution in standard format for pcf2deer 
% from the distribution p0(r0) that may not have equidistant sampling
% points and may not completely overlap with the range of vector r
%
% (c) G. Jeschke, 2003
%

% Standard resolution 2000 points


dr=r(2)-r(1);
rmit=min(r0);
if rmit<min(r), rmit=min(r); end
rmat=max(r0);
if rmat>max(r), rmat=max(r); end

rmip=round(rmit/dr);
rmap=round(rmat/dr);
rax=linspace(rmip*dr,rmap*dr,rmap-rmip+1);
distr0=interp1(r0,p0,rax,'pchip',0);
distr=0*r;
rbas=round(min(r)/dr);
distr(rmip-rbas+1:rmap-rbas+1)=distr0;
