function deer_fit = fit_deer(t,vr,r,distr,Pake_base)

ff = form_factor_deer(r/10,distr,t,Pake_base);
k0 = 0.2/t(end);
par0 = [k0 0.25];
par = fminsearch(@rms_bckg,par0,[],t,vr/max(vr),ff);
deer_fit = sim_bckg(par,t,ff);

function rmsd = rms_bckg(v,texp,vexp,simff)

sim = sim_bckg(v,texp,simff);
rmsd = sqrt(sum((vexp-sim).^2));

function [sim,ff,bckg] = sim_bckg(v,texp,simff)

kdec = v(1);
moddepth = v(2);

bckg = real(exp(-(kdec*texp)));
ff = moddepth*simff + (1-moddepth);
sim = ff.*bckg;
bckg = bckg*(1-moddepth);