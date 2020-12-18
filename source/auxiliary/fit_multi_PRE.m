function fom = fit_multi_PRE(v,predictions,parameters,opt)
% Computes the SAXS and SANS curves for an ensemble and their sum of
% chi^2 with respect to experimental curves
%
% v             vector of populations for the ensemble members
% predictions   array [m,n], where m is the number experimental PREs and
%               n-1 the number of ensemble members, n-1 = length(v)
%               1st column    : experimental PRE data
%               columns 2... n: predicte relaxation enhancements Gamma2
% parameters    parameters for blocks of data, array of struct with fields
%               .td             total INEPT time [s]
%               .R2dia          relaxation rate for diamagnetic sample,
%                               [Hz]
%               .range          first and last index of block into array
%                               predictions
% opt           information for plot/display
%               .interactive    if true, interactive mode, defaults to false
%               .plot_axes      axes object for plotting, no plotting if empty,
%                               defaults to no plotting
%               .update_nr      number of iterations between plot updates, defaults
%                               to 1000
%               .old_size       old ensemble size, defaults to zero
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
n = 0;
for kr = 1:length(parameters)
    td = parameters(kr).td;
    R2dia = parameters(kr).R2dia;
    range = parameters(kr).range;
    n = n + 1 + range(2) - range(1);
    coeff = v/sum(v);
    Gamma2 = predictions(range(1):range(2),2:end)*coeff';
    all_pre = R2dia*exp(-td*Gamma2)./(Gamma2+R2dia);
    fom = fom + sum((all_pre-predictions(range(1):range(2),1)).^2);
end
fom = fom/n;

if isfield(opt,'interactive') && opt.interactive
    if mod(call_count,opt.update_nr) == 0
%         figure(777); clf;
%         hold on;
%         plot(predictions(range(1):range(2),1),'k.');
%         plot(all_pre,'ro');
%         figure(1);
        axis(opt.plot_axes);
        cla; hold on;
        plot(1:opt.old_size,v(1:opt.old_size),'o','Color',[0.6,0,0]);
        plot(opt.old_size+1:length(v),v(opt.old_size+1:length(v)),'.','Color',[0.6,0,0]);
        coeff = v(v > opt.threshold*max(v));
        title(sprintf('PRE mean square dev. (%i): %5.3f with %i conformers',call_count/1000,fom,length(coeff)));
        drawnow
    end
    call_count = call_count + 1;
end