function fom = fit_multi_SAS(v,fit,opt)
% Computes the SAXS and SANS curves for an ensemble and their sum of
% chi^2 with respect to experimental curves
%
% v     vector of populations for the ensemble members
% fit   cell vector of distance distribution information, each cell
%       contains data for one SAS restraint, an array [m,n], where
%       m   number of points in angle domain
%       n-3 number of ensemble members, must match length of v
%       1st column    : angle axis (not used here)
%       2nd column    : SAS curve distance distribution
%       3rd column    : SAS error curve
%       columns 4... n: predicted SAS curves for conformers
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

fom = 0;
for k = 1:length(fit)
    data = fit{k};
    curve = data(:,2);
    errors = data(:,3);
    [m,n] = size(data);
    basis = data(:,4:n);
    coeff = v(1:n-3);
    sim = basis*coeff' + v(n-3+k);
    sc = sum(curve.*sim)/sum(sim.*sim);
    chi2 = sum(((curve-sc*sim)./errors).^2)/(m-1);
    fom = fom + chi2;
end

if isfield(opt,'interactive') && opt.interactive
    if mod(call_count,opt.update_nr) == 0
        axis(opt.plot_axes);
        cla; hold on;
        plot(1:opt.old_size,v(1:opt.old_size),'o','Color',[0.6,0,0]);
        plot(opt.old_size+1:length(v),v(opt.old_size+1:length(v)),'.','Color',[0.6,0,0]);
        coeff = v(v > opt.threshold*max(v));
        title(sprintf('SAS chi^2(%i): %5.3f with %i conformers',call_count/1000,fom,length(coeff)));
        drawnow
    end
    call_count = call_count + 1;
end