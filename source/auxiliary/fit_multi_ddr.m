function fom = fit_multi_ddr(v,fit,opt)
% Computes the distance distributions for an ensemble and their geometric
% mean overlap deficiency with respect to restraint distributions
%
% v     vector of populations for the ensemble members
% fit   cell vector of distance distribution information, each cell
%       contains data for one DEER restraint, an array [m,n], where
%       m   number of points in distance domain
%       n-2 number of ensemble members, must match length of v
%       1st column    : distance axis (not used here)
%       2nd column    : restraint distance distribution
%       columns 3... n: distance distributions for ensemble members
% opt   information for plot/display
%       .interactive    if true, interactive mode, defaults to false
%       .plot_axes      axes object for plotting, no plotting if empty,
%                       defaults to no plotting
%       .update_nr      number of iterations between plot updates, defaults
%                       to 1000
%       .old_size       old ensemble size, defaults to zero
%
% G. Jeschke, 2019-2020

if ~exist('opt','var') || isempty(opt)
    opt.plot_axes = [];
    opt.update_nr = 1000;
    opt.interactive = false;
    opt.old_size = 0;
end

if ~isfield(opt,'plot_axes')
    opt.plot_axes = [];
end

if ~isfield(opt,'update_nr')
    opt.update_nr = 1000;
end

if ~isfield(opt,'interactive')
    opt.interactive = false;
end

if ~isfield(opt,'old_size')
    opt.old_size = 0;
end

persistent call_count
if isempty(call_count)
    call_count = 0;
end

fom = 1;
for k = 1:length(fit)
    data = fit{k};
    restraint = data(:,2);
    [~,n] = size(data);
    basis = data(:,3:n);
    coeff = v;
    sim = basis*coeff';
    sim = sim/sum(sim);
    restraint = restraint/sum(restraint);
    overlap = sum(min([restraint';sim']));
    fom = fom * overlap;
end

if fom < 0
    fom = 0;
end

fom = 1 - fom^(1/length(fit));

if isfield(opt,'interactive') && opt.interactive
    if mod(call_count,opt.update_nr) == 0
        axis(opt.plot_axes);
        cla; hold on;
        plot(1:opt.old_size,v(1:opt.old_size),'o','Color',[0.6,0,0]);
        plot(opt.old_size+1:length(v),v(opt.old_size+1:length(v)),'.','Color',[0.6,0,0]);
        coeff = v(v > opt.threshold*max(v));
        title(sprintf('Overlap deficiency(%i): %5.3f with %i conformers',call_count/1000,fom,length(coeff)));
        drawnow
    end
    call_count = call_count + 1;
end