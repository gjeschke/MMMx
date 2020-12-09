function fom = fit_multi_restraints(v,predictions,normalize,opt)
%
% fom = FIT_MULTI_RESTRAINTS(v,predictions,normalize,opt)
%
% Computes integrative fit loss for combined fitting of populations in an
% ensemble that should fulfill DEER and SAXS/SANS restraints
%
% v             vector of populations for the ensemble members
% predictions   struct of cell vectors of predictions per restraint 
%               and conformer for all restraints in the order:
%               .ddr    distance distribution restraints
%               .sas    small-angle scattering curves
%               .pre    paramagnetic relaxation enhancements, to be
%                       implemented
% opt   information for plot/display
%       .interactive    if true, interactive mode, defaults to false
%       .plot_axes      axes object for plotting, no plotting if empty,
%                       defaults to no plotting
%       .update_nr      number of iterations between plot updates, defaults
%                       to 1000
%       .old_size       old ensemble size, defaults to zero
%
% predictions:
% ddr
%       cell vector of distance distribution information, each cell
%       contains data for one DEER restraint, an array [m,n], where
%       m   number of points in distance domain
%       n-2 number of ensemble members, must match length of v
%       1st column    : distance axis (not used here)
%       2nd column    : restraint distance distribution
%       columns 3... n: distance distributions for ensemble members
% sas
%       cell vector of SAS information, each cell
%       contains data for one SAS restraint, an array [m,n], where
%       m   number of points in angle domain
%       n-3 number of ensemble members, must match length of v
%       1st column    : angle axis (not used here)
%       2nd column    : SAS curve distance distribution
%       3rd column    : SAS error curve
%       columns 4... n: predicted SAS curves for conformers
% pre
%       cell vector of SAS information, each cell
%       contains data for one PRE restraint, an array [m,n], where
%       m   number of assigned protons with PRES
%       n-1 number of ensemble members, must match length of v
%       1st column    : residue numbers (not used here)
%       2nd column    : experimental PREs
%       columns 3... n: predicted PREs for conformers
%
% normalize
%       best-fit figures of merit for individual restraint sets, NaN values
%       indicate that the set does not exist or is not fitted, structure
%       .sas    sum of chi^2 for all SAXS/SANS curves
%       .ddr    geometric mean overlap deficiency for all distance
%               distributions
%       .pre    RMSD of all PREs
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
n_sets = 0;
if ~isnan(normalize.sas)
    fom_sas = fit_multi_SAS(v,predictions.sas);
    fom = fom + fom_sas/normalize.sas;
    n_sets = n_sets + 1;
end
if ~isnan(normalize.ddr)
    [~,n] =size(predictions.ddr{1});
    fom_ddr = fit_multi_ddr(v(1:n-2),predictions.ddr);
    fom = fom + fom_ddr/normalize.ddr;
    n_sets = n_sets + 1;
end
if ~isnan(normalize.pre)
    [~,n] =size(predictions.pre);
    fom_pre = fit_multi_pre(v(1:n-2),predictions.pre);
    fom = fom + fom_pre/normalize.pre;
    n_sets = n_sets + 1;
end

fom = fom/n_sets - 1;

if isfield(opt,'interactive') && opt.interactive
    if mod(call_count,opt.update_nr) == 0
        axis(opt.plot_axes);
        cla; hold on;
        plot(1:opt.old_size,v(1:opt.old_size),'o','Color',[0.6,0,0]);
        plot(opt.old_size+1:length(v),v(opt.old_size+1:length(v)),'.','Color',[0.6,0,0]);
        coeff = v(v > opt.threshold*max(v));
        title(sprintf('Loss of merit(%i): %5.3f with %i conformers',call_count/1000,fom,length(coeff)));
        drawnow
    end
    call_count = call_count + 1;
end