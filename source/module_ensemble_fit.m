function [entity,exceptions,failed,restraints,fit_task] = module_ensemble_fit(control,logfid,entity0)
%
% MODULE_ENSEMBLE_FIT Runs the ensemble fit module
%
%   [exceptions,entity] = MODULE_ENSEMBLE_FIT(control,logfid)
%   Fits an ensemble by integrating distance distribution, small-angle
%   scattering, and paramagnetic relaxation enhancement restraints
%   can also be used to compute restraint fulfilment for an existing
%   ensemble
%
% INPUT
% control       control structure with fields
%               .name           'ensemblefit', unused
%               .options        options struct
%               .directives     array of struct with directives
%               .entity_output  number of the entity to be output
% logfid        file handle for log file
% entity0       optional MMMx entity with rigid-body arrangments, used only
%               with 'expand' directive
%
% OUTPUT
% entity        output entity, empty in NOFIT mode
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
% failed        flag indicated whther the module failed
% restraints    restraint structure with fit values
% fit_task      the fit_task defined by control together with evaluated
%               restraints for individual conformers
%
% SUPPORTED DIRECTIVES
%
% initial       initial ensemble (file name in .options{1})
% nofit         suppresses fitting if present, only restraint evaluation
% interactive   if present, information on progress is displayed
% addpdb        PDB files of conformers to be added, file name that can
%               contain wildcards 
% plot          generate plots
% blocksize     size of a block of conformers fitted simultaneously,
%               defaults to 100
% save          file name for saving ensemble, if not given, the ensemble
%               is saved as ensemble.ens to the current directory
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
exceptions = {[]};
warnings = 0;
if ~exist('entity0','var')
    entity0 = [];
end
entity = [];
failed = false;
expand_rba = false;
fname_basis = 'MMMx_expand';

% set defaults

pre_default_td = 10e-3; % 10 ms default total INEPT time
pre_default_tr = 20; % 20 ns default rotational correlation time
pre_default_taui = 0.5; % 0.5 ns default label internal correlation time
pre_max_Gamma2 = 170;

run_fit = true;
plot_result = false;
initial_ensemble = '';
added_conformers = '';
opt.threshold = 0.01; % conformers wih population below this threshold are discarded

opt.interactive = false; % no interactive plotting during fits
opt.blocksize = 100;
opt.skip_restraints = false;

fit_mean_distances = false;

outname = 'ensemble.ens';

sas_options.lm = 20;
sas_options.fb = 18;
sas_options.delete = true;
sas_options.cst = false;
sas_options.err = true;
sas_options.crysol3 = false; 

restraints.ddr(1).labels{1} = '';
restraints.pre(1).label = '';
restraints.sas(1).data{1} = '';

ddr_poi = 0;
sas_poi = 0;
pre_poi = 0;
% read restraints
for d = 1:length(control.directives)
    switch lower(control.directives(d).name)
        case 'nofit'
            run_fit = false;
        case 'expand'
            expand_rba = true;
            if ~isempty(control.directives(d).options{1})
                fname_basis = control.directives(d).options{1};
            end
        case 'plot'
            plot_result = true;
        case 'initial'
            initial_ensemble = control.directives(d).options{1};
        case 'interactive'
            opt.interactive = true;
        case 'blocksize'
            opt.blocksize = str2double(control.directives(d).options{1});
        case 'addpdb'
            added_conformers = control.directives(d).options{1};
        case 'save'
            outname = control.directives(d).options{1};
        case 'rmean'
            fit_mean_distances = true;
        case 'saxs'
            sas_poi = sas_poi + 1;
            restraints.sas(sas_poi).type = 'saxs';
            restraints.sas(sas_poi).data = control.directives(d).options{1};
            restraints.sas(sas_poi).options = sas_options;
            if length(control.directives(d).options) > 1
                if strcmp(control.directives(d).options{2},'crysol3')
                    restraints.sas(sas_poi).options.crysol3 = true;
                end
            end
        case 'sans'
            sas_poi = sas_poi + 1;
            restraints.sas(sas_poi).type = 'sans';
            restraints.sas(sas_poi).data = control.directives(d).options{1};
            restraints.sas(sas_poi).options = sas_options;
            if length(control.directives(d).options) > 1
                restraints.sas(sas_poi).illres = control.directives(d).options{2};
            end
            if length(control.directives(d).options) > 2
                restraints.sas(sas_poi).D2O = str2double(control.directives(d).options{3});
                restraints.sas(sas_poi).options.D2O = restraints.sas(sas_poi).D2O;
            end
        case 'ddr'
            ddr_poi = ddr_poi + 1; % increase ddr block counter
            fprintf(logfid,'ddr %s',control.directives(d).options{1}); % echo directive to log file
            restraints.ddr(ddr_poi).labels{1} = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % different labels
                restraints.ddr(ddr_poi).labels{2} = control.directives(d).options{2};
                fprintf(logfid,' %s\n',control.directives(d).options{2});
            else % same label at both sites
                restraints.ddr(ddr_poi).labels{2} = control.directives(d).options{1};
                fprintf(logfid,'\n');
            end
            [nr,args] = size(control.directives(d).block);
            if nr > 0 % at least one restraint in this block
                restraints.ddr(ddr_poi).site1{nr} = ' ';
                restraints.ddr(ddr_poi).site2{nr} = ' ';
                restraints.ddr(ddr_poi).r(nr) = 0;
                restraints.ddr(ddr_poi).sigr(nr) = 0;
                restraints.ddr(ddr_poi).file{nr} = '*';
            else % empty block
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_ensemble_fit:empty_ddr_block', 'ddr block %i is empty',d);
                record_exception(exceptions{warnings},logfid);
            end
            if args < 3 % misformed restraint line
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_ensemble_fit:misformed_ddr', 'ddr restraint has less than three arguments');
                record_exception(exceptions{warnings},logfid);
                failed = true;
                return
            end
            for kr = 1:nr
                restraints.ddr(ddr_poi).site1{kr} = control.directives(d).block{kr,1};
                restraints.ddr(ddr_poi).site2{kr} = control.directives(d).block{kr,2};
                restraints.ddr(ddr_poi).file{kr} = '';
                arg3 = control.directives(d).block{kr,3};
                if arg3(1) == '@'
                    restraints.ddr(ddr_poi).file(kr) = arg3(2:end);
                    restraints.ddr(ddr_poi).r(kr) = [];
                    restraints.ddr(ddr_poi).sigr(kr) = [];
                else
                    restraints.ddr(ddr_poi).r(kr) = str2double(arg3);
                    restraints.ddr(ddr_poi).sigr(kr) = str2double(control.directives(d).block{kr,4});
                end
                if args > 4
                    arg5 = control.directives(d).block{kr,5};
                    if arg5(1) == '@'
                        restraints.ddr(ddr_poi).file{kr} = arg5(2:end);
                    end
                end
                for karg = 1:args
                    fprintf(logfid,'  %s',control.directives(d).block{kr,karg});
                end
                fprintf(logfid,'\n');
            end
            fprintf(logfid,'\n\n');
        case {'pre','prerates'}
            pre_poi = pre_poi + 1; % increase pre block counter
            fprintf(logfid,'pre %s',control.directives(d).options{1}); % echo directive to log file
            restraints.pre(pre_poi).label = control.directives(d).options{1};
            if strcmpi(control.directives(d).name,'prerates')
                restraints.pre(pre_poi).fit_rates = true;
            else
                restraints.pre(pre_poi).fit_rates = false;
            end                
            for karg = 2:length(control.directives(d).options) % echo additional arguments
                fprintf(logfid,' %s',control.directives(d).options{karg});
            end
            fprintf(logfid,'\n');
            restraints.pre(pre_poi).site = control.directives(d).options{2};
            restraints.pre(pre_poi).larmor = str2double(control.directives(d).options{3});
            if length(control.directives(d).options) > 3
                restraints.pre(pre_poi).td = 1e-3*str2double(control.directives(d).options{4}); % convert from ms to s
            else
                restraints.pre(pre_poi).td = pre_default_td;
            end
            if length(control.directives(d).options) > 4
                restraints.pre(pre_poi).R2dia = str2double(control.directives(d).options{5}); 
            else
                restraints.pre(pre_poi).R2dia = 0; % there is no default, if R2dia is missing
            end
            if length(control.directives(d).options) > 5
                restraints.pre(pre_poi).taui = 1e-9*str2double(control.directives(d).options{6}); 
            else
                restraints.pre(pre_poi).taui = 1e-9*pre_default_taui; % set default
            end
            if length(control.directives(d).options) > 6
                restraints.pre(pre_poi).taur = str2double(control.directives(d).options{7}); 
            else
                restraints.pre(pre_poi).taur = NaN; % computation will be attempted, if taur is missing
            end
            if length(control.directives(d).options) > 7
                restraints.pre(pre_poi).max_Gamma2 = str2double(control.directives(d).options{8}); 
            else
                restraints.pre(pre_poi).max_Gamma2 = pre_max_Gamma2; % enhancement, where signal is broadened beyond detection
            end
            [nr,args] = size(control.directives(d).block);
            if nr > 0 % at least one restraint in this block
                restraints.ddr(pre_poi).residue{nr} = ' ';
                restraints.ddr(pre_poi).data = ones(2,nr);
                restraints.ddr(pre_poi).data(2,:) = 0;
            else % empty block
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_ensemble_fit:empty_pre_block', 'pre block %i is empty',d);
                record_exception(exceptions{warnings},logfid);
            end
            if args < 2 % misformed restraint line
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_ensemble_fit:misformed_pre', 'pre restraint has less than three arguments');
                record_exception(exceptions{warnings},logfid);
                failed = true;
                return
            end
            for kr = 1:nr
                restraints.pre(pre_poi).residue{kr} = control.directives(d).block{kr,1};
                restraints.pre(pre_poi).data(kr,1) = str2double(control.directives(d).block{kr,2});
                if args > 2
                    restraints.pre(pre_poi).data(kr,2) = str2double(control.directives(d).block{kr,3});
                else
                    restraints.pre(pre_poi).data(kr,2) = 0;
                end                    
                for karg = 1:args
                    fprintf(logfid,'  %s',control.directives(d).block{kr,karg});
                end
                fprintf(logfid,'\n');
            end
            fprintf(logfid,'\n\n');
        otherwise
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_ensemble_fit:unknown_directive',...
                'directive %s is unknown',lower(control.directives(d).name));
            record_exception(exceptions{warnings},logfid);
    end
end

restraints.ddr = restraints.ddr(1:ddr_poi);
restraints.sas = restraints.sas(1:sas_poi);
restraints.pre = restraints.pre(1:pre_poi);

if expand_rba
    if isfield(entity0,'rba_populations')
        C = length(entity0.rba_populations);
    else
        warnings = warnings + 1;
        exceptions{warnings} = MException('module_ensemble_fit:no_rba',...
                'Rigid-body expansion switched off, as there are no rigid bodies in input entity');   
        record_exception(exceptions{warnings},logfid);
        expand_rba = false;
    end
end

[initial_files,pop] = rd_ensemble_definition(initial_ensemble);
[~,~,ext] = fileparts(added_conformers);
if isempty(ext)
    added_conformers = strcat(added_conformers,'.pdb');
end
added_files = dir(added_conformers); % find all files that match the pattern

if ~expand_rba
    % make file list of conformers
    C = length(initial_files) + length(added_files);
    expand_str = '';
else
    expand_str = ' by rigid-body arrangment expansion';
end

fprintf(logfid,'Ensemble fit with %i conformers%s\n',C,expand_str);

C0 = length(initial_files); % number of conformers that must be included

if C == 0
    warnings = warnings + 1;
    exceptions{warnings} = MException('module_ensemble_fit:no_conformers', 'No conformers are specified or found.');
    record_exception(exceptions{warnings},logfid);
    failed = true;
    fit_task = [];
    return
end

fit_task.file_list = cell(1,C);
for c = 1:length(initial_files)
    fit_task.file_list{c} = initial_files(c).name;
end
for c = 1:length(added_files)
    fit_task.file_list{length(initial_files)+c} = added_files(c).name;
end

if C > length(initial_files)
    fit_task.pop = ones(1,C)/C;
else
    fit_task.pop = pop;
end

% evaluate restraints for individual conformers

% distance distribution restraints
options.rmin = 10; % minimal distance
options.rmax = 150; % maximal distance
options.resolution = 0.5;
fit_task.r_axis = options.rmin:options.resolution:options.rmax;

nr = 0;
for ddr_poi = 1:length(restraints.ddr)
    nr = nr + length(restraints.ddr(ddr_poi).site1);
end
if ~isempty(restraints.ddr)
    fit_task.ddr(nr).sim_distr = zeros(1,length(fit_task.r_axis));
    fit_task.ddr(nr).exp_distr = [];
    fit_task.ddr(nr).exp_distr_lb = [];
    fit_task.ddr(nr).exp_distr_ub = [];
    fit_task.ddr(nr).fit_distr = [];
end
kft = 0;
for ddr_poi = 1:length(restraints.ddr)
    for kr = 1:length(restraints.ddr(ddr_poi).site1)
        kft = kft + 1;
        if ~isempty(restraints.ddr(ddr_poi).file{kr})
            exp_data = load(restraints.ddr(ddr_poi).file{kr});
            rax_exp = 10*exp_data(:,1)'; % conversion to Angstroem
            distr_exp = exp_data(:,2)';
            lb_exp = exp_data(:,3)';
            ub_exp = exp_data(:,4)';
            fit_task.ddr(kft).exp_distr = interp1(rax_exp,distr_exp,fit_task.r_axis,'pchip',0);
            sc = 1/sum(fit_task.ddr(kft).exp_distr);
            fit_task.ddr(kft).exp_distr = sc*fit_task.ddr(kft).exp_distr;
            fit_task.ddr(kft).exp_distr_lb = sc*interp1(rax_exp,lb_exp,fit_task.r_axis,'pchip',0);
            fit_task.ddr(kft).exp_distr_ub = sc*interp1(rax_exp,ub_exp,fit_task.r_axis,'pchip',0);
        end
        if ~isempty(restraints.ddr(ddr_poi).r(kr))
            rmean = restraints.ddr(ddr_poi).r(kr);
            sigr = restraints.ddr(ddr_poi).sigr(kr);
            rax = fit_task.r_axis;
            sim_distr = exp(-(rax-rmean).^2/(2*sigr^2));
            fit_task.ddr(kft).sim_distr = sim_distr/sum(sim_distr);
        end
        fit_task.ddr(kft).distr = zeros(C,length(fit_task.r_axis));
    end
end

nr_pre = 0;
for pre_poi = 1:length(restraints.pre)
    nr_pre = nr_pre + length(restraints.pre(pre_poi).residue);
end
if ~isempty(restraints.pre)
    fit_task.pre(nr_pre).exp_data = 1;
    fit_task.pre(nr_pre).fit_data = [];
    fit_task.pre(nr_pre).residue = '';
end
kft = 0;
for pre_poi = 1:length(restraints.pre)
    for kr = 1:length(restraints.pre(pre_poi).residue)
        kft = kft + 1;
        fit_task.pre(kft).exp_data = restraints.pre(pre_poi).data(kr,:); 
        if restraints.pre(pre_poi).fit_rates && fit_task.pre(kft).exp_data(1) > restraints.pre(pre_poi).max_Gamma2
            fit_task.pre(kft).exp_data(1) = restraints.pre(pre_poi).max_Gamma2;
            fprintf(logfid,'Warning: PRE rate for %s@%s.%s truncated to maximum value of %6.1f/s\n',...
                restraints.pre(pre_poi).label,restraints.pre(pre_poi).site,restraints.pre(pre_poi).residue{kr},...
                restraints.pre(pre_poi).max_Gamma2);
        end
        fit_task.pre(kft).residue = restraints.pre(pre_poi).residue{kr};
    end
end

if ~isempty(restraints.sas)
    fit_task.sas(length(restraints.sas)).fits = [];
    sas_initialized = false(1,length(restraints.sas));
end
fit_task.ddr_valid = true(1,nr); % for logical indexing of restraints
fit_task.pre_valid = true(1,nr_pre); % for logical indexing of restraints
if ~isempty(restraints.pre)
    pre_parameters(length(restraints.pre)).td = [];
    pre_parameters(length(restraints.pre)).R2dia = [];
end
for nr_pre = 1:length(restraints.pre)
    pre_parameters(nr_pre).range = [1e6,-1];
end
valid_conformers = true(1,C); % for logical indexing of conformers
for c = 1:C
    if expand_rba
        entity1 = get_rba(entity0,c);
        fname = sprintf('%s_rba_%i.pdb',fname_basis,c);
        put_pdb(entity1,fname);
        fit_task.file_list{c} = fname;
    else
        fname = fit_task.file_list{c};
    end
    fprintf(1,'Working on conformer: %s\n',fname);
    poi = strfind(fname,'.pdb');
    fname0 = fname;
    if ~isempty(poi)
        fname = fname(1:poi-1);
    end
    prop_file = strcat(fname,'.prop.mat');
    if exist(prop_file,'file')
        p = load(prop_file);
        properties = p.properties;
        if isfield(properties,'entity')
            entity = properties.entity;
        else
            entity = get_pdb(fname0);
        end
    else
        clear properties
        properties.initialized = true;
        entity = get_pdb(fname0);
    end
    block_offset = 0;
    for ddr_poi = 1:length(restraints.ddr)
        label1 = restraints.ddr(ddr_poi).labels{1};
        label2 = restraints.ddr(ddr_poi).labels{2};
        for kr = 1:length(restraints.ddr(ddr_poi).site1)
            fit_task.ddr(block_offset+kr).assignment = [ddr_poi,kr];
            site1 = restraints.ddr(ddr_poi).site1{kr};
            site2 = restraints.ddr(ddr_poi).site2{kr};
            is_known = false;
            if isfield(properties,'ddr')
                n_comp = length(properties.ddr);
                for k_comp = 1:n_comp
                    if strcmp(properties.ddr(k_comp).site1,site1) && ...
                            strcmp(properties.ddr(k_comp).site2,site2) && ...
                            strcmp(properties.ddr(k_comp).label1,label1) && ...
                            strcmp(properties.ddr(k_comp).label2,label2)
                        ddr_exceptions = [];
                        distr = properties.ddr(k_comp).distr;
                        is_known = true;
                        break
                    end                    
                end
            else
                n_comp = 0;
            end
            if ~is_known
                [~,distr,entity,ddr_exceptions] = distance_distribution(entity,site1,label1,site2,label2,options);
            end
            if ~isempty(ddr_exceptions) && ~isempty(ddr_exceptions{1})
                for ex = 1:length(ddr_exceptions)
                    warnings = warnings + 1;
                    exceptions{warnings} = ddr_exceptions{ex};
                end
            elseif ~is_known
                n_comp = n_comp + 1;
                properties.ddr(n_comp).site1 = site1;
                properties.ddr(n_comp).site2 = site2;
                properties.ddr(n_comp).label1 = label1;
                properties.ddr(n_comp).label2 = label2;
                properties.ddr(n_comp).distr = distr;
            end
            if ~isempty(distr)
                fit_task.ddr(block_offset+kr).distr(c,:) = distr/sum(distr);
            else
                fprintf(logfid,'Warning: ddr between sites %s and %s cannot be evaluated for conformer %i\n',site1,site2,c);
                if opt.skip_restraints % if requested, skip restraints that cannot be evaluated for all conformers
                    fit_task.ddr_valid(block_offset+kr) = false;
                else % otherwise skip conformers for which not all restraints can be evaluated
                    valid_conformers(c) = false;
                end
            end
        end
        block_offset = block_offset + length(restraints.ddr(ddr_poi).site1);
    end
    block_offset = 0;
    for pre_poi = 1:length(restraints.pre)
        label = restraints.pre(pre_poi).label;
        site = restraints.pre(pre_poi).site;
        site_list(length(restraints.pre(pre_poi).residue)).chain = ''; %#ok<AGROW>
        site_list(length(restraints.pre(pre_poi).residue)).residue = []; %#ok<AGROW>
        if ~isnan(restraints.pre(pre_poi).taur)
            taur = 1e-9*restraints.pre(pre_poi).taur;
        else
            if isfield(properties,'taur') && ~isempty(properties.taur)
                taur = properties.taur;
            else
                options.N = 1000;
                options.HYDROPRO = true;
                [taur,info] = rotational_correlation_time(entity,'(*)',options);
                if ~info.HYDROPRO
                    taur = 1e-9*pre_default_tr;
                    fprintf(logfid,'Warning: pre %i assumes default rotational correlation time %4.1f ns\n',pre_default_tr);
                end
                properties.taur = taur;
            end
        end
        taui = restraints.pre(pre_poi).taui;
        td = restraints.pre(pre_poi).td;
        R2dia = restraints.pre(pre_poi).R2dia;
        larmor = restraints.pre(pre_poi).larmor;
        pre_parameters(pre_poi).td = td;
        pre_parameters(pre_poi).R2dia = R2dia;
        pre_parameters(pre_poi).fit_rates = restraints.pre(pre_poi).fit_rates;
        pre_parameters(pre_poi).max_Gamma2 = restraints.pre(pre_poi).max_Gamma2;
        for kr = 1:length(restraints.pre(pre_poi).residue)
            fit_task.pre(block_offset+kr).assignment = [pre_poi,kr];
            residue = restraints.pre(pre_poi).residue{kr};
            si = strfind(residue,'(');
            ei = strfind(residue,')');
            if ~isempty(si) && ~isempty(ei)
                chain_str = residue(si+1:ei-1);
            else
                ei = 0;
            end
            site_list(kr).chain = chain_str;
            site_list(kr).residue = str2double(residue(ei+1:end));
        end
        [pre_list,pre_exceptions,entity] = pre_of_entity(entity,site_list,label,site,td,taur,taui,R2dia,larmor);
        if ~isempty(pre_exceptions) && ~isempty(pre_exceptions{1})
            for ex = 1:length(pre_exceptions)
                warnings = warnings + 1;
                exceptions{warnings} = pre_exceptions{ex};
            end
        end
        for kr = 1:length(site_list)
            found = false;
            for kcomp = 1:length(pre_list)
                if strcmp(site_list(kr).chain,pre_list(kcomp).chain) && ...
                        site_list(kr).residue == pre_list(kr).residue
                    fit_task.pre(block_offset+kr).Gamma2(c) = pre_list(kr).Gamma2;
                    if ~isnan(pre_list(kr).Gamma2)
                        found = true;
                        break
                    end
                end
            end
            if ~found
                fprintf(logfid,'Warning: pre at sites (%s)%i and %i cannot be evaluated for conformer %i\n',site_list(kr).chain,site_list(kr).residue,pre_list(kr).residue,c);
                if opt.skip_restraints % if requested, skip restraints that cannot be evaluated for all conformers
                    fit_task.pre_valid(block_offset+kr) = false;
                else % otherwise skip conformers for which not all restraints can be evaluated
                    valid_conformers(c) = false;
                end
            end
        end
        block_offset = block_offset + length(restraints.pre(pre_poi).residue);
    end
    fit_task.sas_valid = false(1,length(restraints.sas)); % for logical indexing of restraints
    for sas_poi = 1:length(restraints.sas)
        is_known = false;
        if isfield(properties,'sas')
            n_comp = length(properties.sas);
            for k_comp = 1:n_comp
                if strcmp(restraints.sas(sas_poi).type,properties.sas(k_comp).type) && ...
                        strcmp(restraints.sas(sas_poi).data,properties.sas(k_comp).data)
                    fit = properties.sas(k_comp).fit;
                    chi2 = properties.sas(k_comp).chi2;
                    is_known = true;
                    break
                end
            end
            if ~is_known
                datafile = restraints.sas(sas_poi).data;
                options = restraints.sas(sas_poi).options;
                switch restraints.sas(sas_poi).type
                    case 'saxs'
                        [fit,chi2] = fit_SAXS(datafile,fit_task.file_list{c},options);
                        figure(3); clf; hold on;
                        plot(fit(:,1),fit(:,2),'k');
                        plot(fit(:,1),fit(:,3),'r');
                        title(sprintf('chi^2 = %6.3f\n',chi2));
                        fprintf(1,'SAXS chi^2 = %6.3f\n',chi2);
                        drawnow
                    case 'sans'
                        illres = restraints.sas(sas_poi).illres;
                        [fit,chi2] = fit_SANS(datafile,fit_task.file_list{c},illres,options);
                end
            end
            if isempty(fit)
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_ensemble_fit:sas_curve_fitting_failed', 'SAS curve %s could not be fitted for conformer %s',datafile,fit_task.file_list{c});
                fprintf(logfid,'Warning: sas %i cannot be evaluated for conformer %i\n',sas_poi,c);
                valid_conformers(c) = false;
            else
                if ~sas_initialized(sas_poi)
                    [m,~] = size(fit);
                    fit_task.sas(sas_poi).fits = zeros(m,3+C);
                    fit_task.sas(sas_poi).fits(:,1:2) = fit(:,1:2);
                    fit_task.sas(sas_poi).fits(:,3) = fit(:,4);
                    fit_task.sas(sas_poi).fits(:,4) = fit(:,3);
                    sas_initialized(sas_poi) = true;
                else
                    fit_task.sas(sas_poi).fits(:,3+c) = fit(:,3);
                end
                fit_task.sas_valid(sas_poi) = true; % if at least one conformer can be computed
            end
        else
            n_comp = 0;
            datafile = restraints.sas(sas_poi).data;
            options = restraints.sas(sas_poi).options;
            switch restraints.sas(sas_poi).type
                case 'saxs'
                    [fit,chi2] = fit_SAXS(datafile,fit_task.file_list{c},options);
                case 'sans'
                    illres = restraints.sas(sas_poi).illres;
                    [fit,chi2] = fit_SANS(datafile,fit_task.file_list{c},illres,options);
            end
            if isempty(fit)
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_ensemble_fit:sas_curve_fitting_failed', 'SAS curve %s could not be fitted for conformer %s',datafile,fit_task.file_list{c});
                fprintf(logfid,'Warning: sas %i cannot be evaluated for conformer %i\n',sas_poi,c);
                valid_conformers(c) = false;
            else
                if ~sas_initialized(sas_poi)
                    [m,~] = size(fit);
                    fit_task.sas(sas_poi).fits = zeros(m,3+C);
                    fit_task.sas(sas_poi).fits(:,1:2) = fit(:,1:2);
                    fit_task.sas(sas_poi).fits(:,3) = fit(:,4);
                    fit_task.sas(sas_poi).fits(:,4) = fit(:,3);
                    sas_initialized(sas_poi) = true;
                else
                    fit_task.sas(sas_poi).fits(:,3+c) = fit(:,3);
                end
                fit_task.sas_valid(sas_poi) = true; % if at least one conformer can be computed
            end
        end
        if ~is_known
            n_comp = n_comp + 1;
            properties.sas(n_comp).fit = fit;
            properties.sas(n_comp).chi2 = chi2;
            properties.sas(n_comp).data = restraints.sas(sas_poi).data;
            properties.sas(n_comp).type = restraints.sas(sas_poi).type;
        end
    end
    properties.entity = entity;
    save(prop_file,'properties');
end

% restrict to the conformers for which all valid restraints could be
% evaluated
fit_task.file_list = fit_task.file_list(valid_conformers);
fit_task.pop = fit_task.pop(valid_conformers);
fit_task.pop = fit_task.pop/sum(fit_task.pop);
if isfield(fit_task,'sas')
    for kr = 1:length(fit_task.sas)
        fit_task.sas(kr).all_fits = fit_task.sas(kr).fits(:,4:end);
        sas_fits = fit_task.sas(kr).fits(:,4:end);
        fit_task.sas(kr).fits = [fit_task.sas(kr).fits(:,1:3) sas_fits(:,valid_conformers)];
    end
end
if isfield(fit_task,'pre')
    for kr = 1:length(fit_task.pre)
        fit_task.pre(kr).Gamma2 = fit_task.pre(kr).Gamma2(valid_conformers);
    end
end
C = length(fit_task.pop);

fprintf(logfid,'\n');

fit_task.remaining_conformers = 1:C;
fit_task.ensemble_populations = fit_task.pop;

% prepare DEER fit by arranging fit target and distance distributions
% predicted for all conformers into matrices for all restraints
nr_ddr = 0;
all_ddr_predictions = cell(1,length(fit_task.ddr_valid));
for kr = 1:length(fit_task.ddr_valid) % all restraints
    if fit_task.ddr_valid(kr) % this is a valid restraint
        nr_ddr = nr_ddr + 1;
        data = zeros(length(fit_task.r_axis),C+2);
        data(:,1) = fit_task.r_axis.';
        fit_task.ddr(kr).distr = fit_task.ddr(kr).distr(valid_conformers,:);
        if ~isempty(fit_task.ddr(kr).exp_distr)
            data(:,2) = fit_task.ddr(kr).exp_distr.';
        else
            data(:,2) = fit_task.ddr(kr).sim_distr.';
        end
        data(:,3:end) = fit_task.ddr(kr).distr.';
        all_ddr_predictions{nr_ddr} = data;
    end
end
all_ddr_predictions = all_ddr_predictions(1:nr_ddr);

% prepare SAS fit by arranging fit target and distance distributions
% predicted for all conformers into matrices for all restraints
nr_sas = 0;
all_sas_predictions = cell(1,length(fit_task.sas_valid));
for kr = 1:length(fit_task.sas_valid) % all restraints
    if fit_task.sas_valid(kr) % this is a valid restraint
        nr_sas = nr_sas + 1;
        all_sas_predictions{kr} = fit_task.sas(kr).fits;
    end
end
all_sas_predictions = all_sas_predictions(1:nr_sas);

% prepare PRE fit by arranging fit target and computed Gamma2 values
% predicted for all conformers into matrices for all restraints
nr_pre = 0;
all_pre_predictions = zeros(length(fit_task.pre_valid),C+2);
for kr = 1:length(fit_task.pre_valid) % all restraints
    if fit_task.pre_valid(kr) % this is a valid restraint
        nr_pre = nr_pre + 1;
        block = fit_task.pre(kr).assignment(1);
        if nr_pre < pre_parameters(block).range(1)
            pre_parameters(block).range(1) = nr_pre;
        end
        if nr_pre > pre_parameters(block).range(2)
            pre_parameters(block).range(2) = nr_pre;
        end
        all_pre_predictions(nr_pre,1:2) = fit_task.pre(kr).exp_data;
        all_pre_predictions(nr_pre,3:C+2) = fit_task.pre(kr).Gamma2;
    end
end
all_pre_predictions = all_pre_predictions(1:nr_pre,:);

loss_of_merit = [];

if run_fit
    % prepare display figure if needed and store axes handle
    if opt.interactive
        figure
        opt.plot_axes = gca;
        opt.old_size = length(initial_files);
    end

    % the following is an iterative run with block size limit
    
    processed = C0; % pointer for processed files, initial ensemble is automatically included
    included = 1:C0; % indices for conformers of initial ensemble
    first_run = true; % even if only the initial ensemble is given, it should be refitted
    
    while processed < C || first_run
    
        first_run = false;
        
        if C0 > opt.blocksize - 10 % increase blocksize if necessary and raise warning 
                                   % process at least 10 additional conformers
            blocksize = C0 + round(opt.blocksize/2);
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_ensemble_fit:blocksize_increased', 'Block size was increased to %i.',blocksize);
            record_exception(exceptions{warnings},logfid);
        else
            blocksize = opt.blocksize;
        end

        % make the list of conformer indices for this block and extract fit task 
        
        conformers = zeros(1,blocksize);
        conformers(1:C0) = included;
        % determine index of last conformer in this block
        last_conformer = processed+blocksize-C0;
        if last_conformer > C
            last_conformer = C;
        end
        % add the new conformers for this block
        for c = processed+1:last_conformer
            conformers(C0+c-processed) = c;
        end
        conformers = conformers(1:C0+last_conformer-processed);
        processed = last_conformer;
        curr_blocksize = length(conformers); % the actual blocksize, as we might have run out of conformers
        
        % initialize figure-of-merit structure with NaNs
        
        fom.ddr = NaN;
        fom.sas = NaN;
        fom.pre = NaN;
        
        nr_sets = 0;
        
        if nr_ddr > 0 % distance distribution fit, if there are such restraints
            clear fit_multi_ddr % initialize iteration counter
            
            curr_ddr_predictions = all_ddr_predictions;
            for kr = 1:length(all_ddr_predictions) % loop over all restraints
                data = all_ddr_predictions{kr};
                data = [data(:,1:2) data(:,2+conformers)];
                curr_ddr_predictions{kr} = data;
            end
            
            if fit_mean_distances % fit only mean distances
                fanonym_ddr = @(v_opt)fit_multi_rmean(v_opt,curr_ddr_predictions,opt);
            else % this is the usual mode
                fanonym_ddr = @(v_opt)fit_multi_ddr(v_opt,curr_ddr_predictions,opt);
            end
            
            l_ddr = zeros(1,curr_blocksize); % lower bound, populations are non-negative
            u_ddr = ones(1,curr_blocksize); % upper bound, populations cannot exceed 1
            v0 = ones(1,curr_blocksize)/curr_blocksize; % start from uniform populations
            
            fit_options = optimoptions('patternsearch',...
                'MaxFunctionEvaluations',10000*length(v0),'MaxIterations',500*length(v0),...
                'StepTolerance',opt.threshold/10);
            
            tic,
            [v,fom_ddr,exitflag,fit_output] = patternsearch(fanonym_ddr,v0,[],[],[],[],l_ddr,u_ddr,[],fit_options);
            ddr_time = toc;
            th = floor(ddr_time/3600);
            ddr_time = ddr_time - 3600*th;
            tmin = floor(ddr_time/60);
            ts = round(ddr_time - 60*tmin);
            
            % store information for multi-restraint fit
            fom.ddr = fom_ddr;
            predictions.ddr = curr_ddr_predictions;
            nr_sets = nr_sets + 1;
            
            fprintf(logfid,'DEER restraints fit took %i h %i min %i s\n',th,tmin,ts);
            switch exitflag
                case 0
                    fprintf(logfid,'Warning: Maximum number of function evaluations or iterations reached. Not converged.\n');
                case 1
                    fprintf(logfid,'Convergence by mesh size.\n');
                case 2
                    fprintf(logfid,'Convergence by change in population.\n');
                case 3
                    if fit_mean_distances
                        fprintf(logfid,'Convergence by mean distance RMSD precision.\n');
                    else
                        fprintf(logfid,'Convergence by ddr overlap deficiency precision.\n');
                    end
                case 4
                    fprintf(logfid,'Convergence by machine precision.\n');
                case -1
                    fprintf(logfid,'ERROR: Optimization terminated by output or plot function.\n');
                case -2
                    fprintf(logfid,'ERROR: No feasible solution found.\n');
                case -3
                    fprintf(logfid,'Warning: No fitting was requested in restraint file.\n');
            end
            fprintf(logfid,'%i iterations and %i function evaluations were performed.\n',fit_output.iterations,fit_output.funccount);
            fprintf(logfid,'Mesh size is at %12.4g, maximum constraint violation at %12.4g.\n',fit_output.meshsize,fit_output.maxconstraint);
            
            coeff_ddr = v/max(v); % populations normalized to their maximum
            above_threshold_ddr = (coeff_ddr >= opt.threshold); % indicies of conformers above the population threshold
            included_ddr = conformers(above_threshold_ddr); % these conformers are kept
            C0_ddr = length(included_ddr);
            if fit_mean_distances
                fprintf(logfid,'Final mean distance rmsd is: %6.3f A with %i conformers\n\n',fom_ddr,C0_ddr);
            else
                fprintf(logfid,'Final ddr overlap deficiency is: %6.3f with %i conformers\n\n',fom_ddr,C0_ddr);
            end
        end
        
        curr_sas_predictions = {};
        predictions.sas = curr_sas_predictions;
        
        if nr_sas > 0 % small-angle scattering curve fit, if there are such restraints
            clear fit_multi_SAS % initialize iteration counter
            
            curr_sas_predictions = all_sas_predictions;
            for kr = 1:length(all_sas_predictions) % loop over all restraints
                data = all_sas_predictions{kr};
                data = [data(:,1:3) data(:,3+conformers)];
                curr_sas_predictions{kr} = data;
            end
            
            fanonym_sas = @(v_opt)fit_multi_SAS(v_opt,curr_sas_predictions,opt);
            
            l_sas = [zeros(1,curr_blocksize)  -0.01*ones(1,length(curr_sas_predictions))]; % lower bound, populations are non-negative, add baseline corrections
            u_sas = [ones(1,curr_blocksize) 0.01*ones(1,length(curr_sas_predictions))]; % upper bound, populations cannot exceed 1, add baseline corrections
            v0 = [ones(1,curr_blocksize)/curr_blocksize zeros(1,length(curr_sas_predictions))]; % start from uniform populations
            
            fit_options = optimoptions('patternsearch',...
                'MaxFunctionEvaluations',10000*length(v0),'MaxIterations',500*length(v0),...
                'StepTolerance',opt.threshold/10);
            tic,
            [v,fom_sas,exitflag,fit_output] = patternsearch(fanonym_sas,v0,[],[],[],[],l_sas,u_sas,[],fit_options);
            sas_time = toc;
            th = floor(sas_time/3600);
            sas_time = sas_time - 3600*th;
            tmin = floor(sas_time/60);
            ts = round(sas_time - 60*tmin);
            
            % store information for multi-restraint fit
            fom.sas = fom_sas;
            predictions.sas = curr_sas_predictions;
            nr_sets = nr_sets + 1;
            
            fprintf(logfid,'SAS curve fit took %i h %i min %i s\n',th,tmin,ts);
            switch exitflag
                case 0
                    fprintf(logfid,'Warning: Maximum number of function evaluations or iterations reached. Not converged.\n');
                case 1
                    fprintf(logfid,'Convergence by mesh size.\n');
                case 2
                    fprintf(logfid,'Convergence by change in population.\n');
                case 3
                    fprintf(logfid,'Convergence by SAS chi^2 precision.\n');
                case 4
                    fprintf(logfid,'Convergence by machine precision.\n');
                case -1
                    fprintf(logfid,'ERROR: Optimization terminated by output or plot function.\n');
                case -2
                    fprintf(logfid,'ERROR: No feasible solution found.\n');
                case -3
                    fprintf(logfid,'Warning: No fitting was requested in restraint file.\n');
            end
            fprintf(logfid,'%i iterations and %i function evaluations were performed.\n',fit_output.iterations,fit_output.funccount);
            fprintf(logfid,'Mesh size is at %12.4g, maximum constraint violation at %12.4g.\n',fit_output.meshsize,fit_output.maxconstraint);
            
            v = v(1:curr_blocksize);
            coeff_sas = v/max(v); % populations normalized to their maximum
            above_threshold_sas = (coeff_sas >= opt.threshold); % indices of conformers above the population threshold
            included_sas = conformers(above_threshold_sas); % these conformers are kept
            C0_sas = length(included_sas);
            fprintf(logfid,'Final SAS chi^2 is: %6.3f with %i conformers\n\n',fom_sas,C0_sas);
        end
        
        if nr_pre > 0 % PRE fit, if there are such restraints
            clear fit_multi_PRE % initialize iteration counter
            
            curr_pre_predictions = [all_pre_predictions(:,1:2) all_pre_predictions(:,2+conformers)];
            
            fanonym_pre = @(v_opt)fit_multi_PRE(v_opt,curr_pre_predictions,pre_parameters,opt);
            
            l_pre = zeros(1,curr_blocksize); % lower bound, populations are non-negative
            u_pre = ones(1,curr_blocksize); % upper bound, populations cannot exceed 1
            v0 = ones(1,curr_blocksize)/curr_blocksize; % start from uniform populations
            
            fit_options = optimoptions('patternsearch',...
                'MaxFunctionEvaluations',10000*length(v0),'MaxIterations',500*length(v0),...
                'StepTolerance',opt.threshold/10);
            tic,
            [v,fom_pre,exitflag,fit_output] = patternsearch(fanonym_pre,v0,[],[],[],[],l_pre,u_pre,[],fit_options);
            pre_time = toc;
            th = floor(pre_time/3600);
            pre_time = pre_time - 3600*th;
            tmin = floor(pre_time/60);
            ts = round(pre_time - 60*tmin);
            
            % store information for multi-restraint fit
            fom.pre = fom_pre;
            predictions.pre.predictions = curr_pre_predictions;
            predictions.pre.parameters = pre_parameters;
            nr_sets = nr_sets + 1;
            
            fprintf(logfid,'PRE fit took %i h %i min %i s\n',th,tmin,ts);
            switch exitflag
                case 0
                    fprintf(logfid,'Warning: Maximum number of function evaluations or iterations reached. Not converged.\n');
                case 1
                    fprintf(logfid,'Convergence by mesh size.\n');
                case 2
                    fprintf(logfid,'Convergence by change in population.\n');
                case 3
                    fprintf(logfid,'Convergence by PRE MSD precision.\n');
                case 4
                    fprintf(logfid,'Convergence by machine precision.\n');
                case -1
                    fprintf(logfid,'ERROR: Optimization terminated by output or plot function.\n');
                case -2
                    fprintf(logfid,'ERROR: No feasible solution found.\n');
                case -3
                    fprintf(logfid,'Warning: No fitting was requested in restraint file.\n');
            end
            fprintf(logfid,'%i iterations and %i function evaluations were performed.\n',fit_output.iterations,fit_output.funccount);
            fprintf(logfid,'Mesh size is at %12.4g, maximum constraint violation at %12.4g.\n',fit_output.meshsize,fit_output.maxconstraint);
            
            coeff_pre = v/max(v); % populations normalized to their maximum
            above_threshold_pre = (coeff_pre >= opt.threshold); % indicies of conformers above the population threshold
            included_pre = conformers(above_threshold_pre); % these conformers are kept
            C0_pre = length(included_pre);
            fprintf(logfid,'Final PRE mean square deviation is: %6.3f with %i conformers\n\n',fom_pre,C0_pre);
        end
        
        loss_of_merit = [];
        if nr_sets > 1 % multi-restraint fit, if there is more than on restraint set
            clear fit_multi_restraints % initialize iteration counter
            
            if fit_mean_distances
                opt.rmean = true;
            else
                opt.rmean = false;
            end
            fanonym_integrative = @(v_opt)fit_multi_restraints(v_opt,predictions,fom,opt);
            
            l = [zeros(1,curr_blocksize)  -0.01*ones(1,length(curr_sas_predictions))]; % lower bound, populations are non-negative, add baseline corrections
            u = [ones(1,curr_blocksize) 0.01*ones(1,length(curr_sas_predictions))]; % upper bound, populations cannot exceed 1, add baseline corrections
            v0 = [ones(1,curr_blocksize)/curr_blocksize zeros(1,length(curr_sas_predictions))]; % start from uniform populations
            
            fit_options = optimoptions('patternsearch',...
                'MaxFunctionEvaluations',50000*length(v0),'MaxIterations',2500*length(v0),...
                'StepTolerance',opt.threshold/10);
            
            tic,
            [v,loss_of_merit,exitflag,fit_output] = patternsearch(fanonym_integrative,v0,[],[],[],[],l,u,[],fit_options);
            integrative_time = toc;
            
            if ~isnan(fom.sas)
                fom_sas = sim_multi_sas(v,curr_sas_predictions);
            end
            v = v(1:end-length(curr_sas_predictions));
            
            if ~isnan(fom.ddr)
                fom_ddr = sim_multi_ddr(v,curr_ddr_predictions);
            end
            
            if ~isnan(fom.pre)
                [fom_pre,pre_fit_type] = sim_multi_pre(v,curr_pre_predictions,pre_parameters);
            end
            
            th = floor(integrative_time/3600);
            integrative_time = integrative_time - 3600*th;
            tmin = floor(integrative_time/60);
            ts = round(integrative_time - 60*tmin);
            
            fprintf(logfid,'Iterative fit took %i h %i min %i s\n',th,tmin,ts);
            switch exitflag
                case 0
                    fprintf(logfid,'Warning: Maximum number of function evaluations or iterations reached. Not converged.\n');
                case 1
                    fprintf(logfid,'Convergence by mesh size.\n');
                case 2
                    fprintf(logfid,'Convergence by change in population.\n');
                case 3
                    fprintf(logfid,'Convergence by loss of merit precision.\n');
                case 4
                    fprintf(logfid,'Convergence by machine precision.\n');
                case -1
                    fprintf(logfid,'ERROR: Optimization terminated by output or plot function.\n');
                case -2
                    fprintf(logfid,'ERROR: No feasible solution found.\n');
                case -3
                    fprintf(logfid,'Warning: No fitting was requested in restraint file.\n');
            end
            fprintf(logfid,'%i iterations and %i function evaluations were performed.\n',fit_output.iterations,fit_output.funccount);
            fprintf(logfid,'Mesh size is at %12.4g, maximum constraint violation at %12.4g.\n',fit_output.meshsize,fit_output.maxconstraint);
            
            coeff_all = v/max(v); % populations normalized to their maximum
            above_threshold_all = (coeff_all >= opt.threshold); % indicies of conformers above the population threshold
            included_all = conformers(above_threshold_all); % these conformers are kept
            C0_all = length(included_all);
            if ~isnan(fom.ddr)
                if fit_mean_distances
                    fprintf(logfid,'Integrative mean distance RMSD is: %6.3f A with %i conformers\n\n',fom_ddr,C0_all);
                else
                    fprintf(logfid,'Integrative ddr overlap deficiency is: %6.3f with %i conformers\n\n',fom_ddr,C0_all);
                end
            end
            if ~isnan(fom.sas)
                fprintf(logfid,'Integrative SAS chi^2 is: %6.3f with %i conformers\n\n',fom_sas,C0_all);            
            end
            if ~isnan(fom.pre)
                fprintf(logfid,'Integrative PRE %s deviation is: %6.4f with %i conformers\n\n',pre_fit_type,fom_pre,C0_all);            
            end
            fprintf(logfid,'Loss of merit is: %6.3f with %i conformers\n\n',loss_of_merit,C0_all);            
        end
        coeff = v/max(v); % populations normalized to their maximum
        above_threshold = (coeff >= opt.threshold); % indicies of conformers above the population threshold
        included = conformers(above_threshold); % these conformers are kept
        fit_task.remaining_conformers = included; % store numbers of remaining conformers
        coeff(~above_threshold) = 0; % discard low-probability conformers
        coeff = coeff/sum(coeff);
        fit_task.ensemble_populations = coeff(above_threshold); % store corresponding populations
        fit_task.pop = coeff; % assign 
        C0 = length(included);
        opt.old_size = C0;
    end
end

fprintf(logfid,'%i conformers remain included\n',C0);          

% save ensemble
ofid = fopen(outname,'wt');
if ofid == -1
    warnings = warnings + 1;
    exceptions{warnings} = MException('module_ensemble_fit:ensemble_not_saved', 'Output file %s could not be written.',outname);
    record_exception(exceptions{warnings},logfid);
else
    fprintf(ofid,'%% %s, %s\n',outname,datetime);
    for c = 1:length(fit_task.remaining_conformers)
        fprintf(ofid,'%s%10.6f\n',fit_task.file_list{fit_task.remaining_conformers(c)},fit_task.ensemble_populations(c));
    end
    fclose(ofid);
end

for kr = 1:nr
    fit_distr = fit_task.ensemble_populations*fit_task.ddr(kr).distr(fit_task.remaining_conformers,:);
    fit_task.ddr(kr).fit_distr = fit_distr;
    restraints.ddr(fit_task.ddr(kr).assignment(1)).fit_distr{fit_task.ddr(kr).assignment(2)} = fit_distr;
end

for kr = 1:nr_sas
%     fitted_curve = fit_task.ensemble_populations*fit_task.sas(kr).all_fits(:,fit_task.remaining_conformers).';
%     fitted_curve = fitted_curve.';
    data = all_sas_predictions{kr};
    data = data(:,[1:3 fit_task.remaining_conformers+3]);
    bsl = 0;
    %    [bsl,chi2] = fminsearch(@sas_baseline,bsl,[],data(:,2),fitted_curve,data(:,3));
    [bsl,chi2] = fminsearch(@sas_baseline,bsl,[],fit_task.ensemble_populations,data);
    curve = data(:,2);
    [~,n] = size(data);
    basis = data(:,4:n);
    sim = basis*fit_task.ensemble_populations' + bsl;
    sc = sum(curve.*sim)/sum(sim.*sim);
%     fitted_curve = fitted_curve + bsl;
%     sc = sum(data(:,2).*fitted_curve)/sum(fitted_curve.^2);            
    fit_task.sas(kr).fitted_curve = sc*sim;
    fit_task.sas(kr).chi2 = chi2;
end

if nr_pre > 0
    kft = 0;
    fom_pre = 0;
    for kr = 1:length(pre_parameters)
        kft0 = 0;
        td = pre_parameters(kr).td;
        R2dia = pre_parameters(kr).R2dia;
        range = pre_parameters(kr).range;
        pre_std = all_pre_predictions(:,2);
        if sum(pre_std == 0) > 0
            do_chi2 = false;
            pre_fit_type = 'mean-square';
        else
            do_chi2 = true;
            pre_fit_type = 'chi-square';
        end
        for nr_pre = range(1):range(2)
            kft = kft + 1;
            kft0 = kft0+1;
            if do_chi2
                pre_std = fit_task.pre(kft).exp_data(2);
            else
                pre_std = 1;
            end
            Gamma2 = sum(fit_task.ensemble_populations.*all_pre_predictions(kft,2+fit_task.remaining_conformers));
            if pre_parameters(kr).fit_rates
                fit_task.pre(kft).fit_data = Gamma2;
            else
                sim_pre = R2dia*exp(-td*Gamma2)/(Gamma2+R2dia);
                fit_task.pre(kft).fit_data = sim_pre;
            end
            fom_pre = fom_pre + ((fit_task.pre(kft).fit_data - fit_task.pre(kft).exp_data(1))/pre_std)^2;
        end
    end
    fom_pre = fom_pre/kft;
end

if isempty(loss_of_merit) % for non-integrative fits, there is no loss of merit
    loss_of_merit = 0;
else
    fprintf(logfid,'\n Final loss of merit: %6.3f\n',loss_of_merit);
    fprintf(1,'\n Final loss of merit: %6.3f\n',loss_of_merit);
end
fit_task.loss_of_merit = loss_of_merit;

% save fit task and restraints

save([outname '.mat'],'restraints','fit_task','exceptions');

if plot_result
    dr = fit_task.r_axis(2) - fit_task.r_axis(1);
    overlap = 1;
    for kr = 1:nr
        if fit_task.ddr_valid(kr)
            figure(kr); clf; hold on
            overlap_G = [];
            overlap_E = [];
            if ~isempty(fit_task.ddr(kr).exp_distr)
                overlap_E = sum(min([fit_task.ddr(kr).fit_distr;fit_task.ddr(kr).exp_distr]));
                if ~isempty(fit_task.ddr(kr).exp_distr_ub)
                    distr_ub = fit_task.ddr(kr).exp_distr_ub;
                    distr_lb = fit_task.ddr(kr).exp_distr_lb;
                    fill([fit_task.r_axis, fliplr(fit_task.r_axis)],dr*[distr_ub, fliplr(distr_lb)],0.75*[1,1,1],'LineStyle','none');
                end
                plot(fit_task.r_axis,dr*fit_task.ddr(kr).exp_distr,'Color',[0.2,0.2,0.2]);
            end
            if ~isempty(fit_task.ddr(kr).sim_distr) && isempty(fit_task.ddr(kr).exp_distr)
                plot(fit_task.r_axis,dr*fit_task.ddr(kr).sim_distr,'Color',[0,0.6,0]);
                overlap_G = sum(min([fit_task.ddr(kr).fit_distr;fit_task.ddr(kr).sim_distr]));
            end
            plot(fit_task.r_axis,dr*fit_task.ddr(kr).fit_distr,'Color',[0.6,0,0]);
            if ~isempty(overlap_E)
                overlap = overlap*overlap_E;
            else
                overlap = overlap*overlap_G;
            end
            site1 = restraints.ddr(fit_task.ddr(kr).assignment(1)).site1{fit_task.ddr(kr).assignment(2)};
            site2 = restraints.ddr(fit_task.ddr(kr).assignment(1)).site2{fit_task.ddr(kr).assignment(2)};
            title_str = sprintf('%s-%s Overlaps:',site1,site2);
            if ~isempty(overlap_E)
                title_str = sprintf('%s exp. %6.3f',title_str,overlap_E);
            elseif ~isempty(overlap_G)
                title_str = sprintf('%s Gauss %6.3f',title_str,overlap_G);
            else
                title_str = sprintf('%s unknown',title_string);
            end
            title(title_str);
            xlabel('Distance (Angstroem)');
            ylabel('Probability density');
        end
    end
    overlap = overlap^(1/nr);
    fprintf(1,'Overlap deficiency: %6.3f\n',1-overlap);
    for kr = 1:nr_sas
        if fit_task.sas_valid(kr)
            figure(kr+1000); clf; hold on
            data = all_sas_predictions{kr};
            fitted_curve = fit_task.sas(kr).fitted_curve;
            plot(data(:,1),data(:,2),'.','MarkerSize',14,'Color',[0.2,0.2,0.2]);
            plot(data(:,1),fitted_curve,'-','LineWidth',2.5,'Color',[0.7,0,0]);
            title(sprintf('chi^2 = %6.3f',fit_task.sas(kr).chi2));
            xlabel('s [Angstroem^{-1}]');
            ylabel('I(s)');
            figure(kr+2000); clf; hold on
            data = all_sas_predictions{kr};
            plot(data(:,1),real(log(data(:,2))),'.','MarkerSize',14,'Color',[0.2,0.2,0.2]);
            plot(data(:,1),real(log(fitted_curve)),'-','LineWidth',2.5,'Color',[0.7,0,0]);
            title(sprintf('chi^2 = %6.3f',fit_task.sas(kr).chi2));
            xlabel('s [Angstroem^{-1}]');
            ylabel('log(I(s))');
            figure(kr+3000); clf; hold on
            data = all_sas_predictions{kr};
            plot(data(:,1),data(:,2)-fitted_curve,'.','MarkerSize',14,'Color',[0,0,0.8]);
            title(sprintf('Fit residual: chi^2 = %6.3f',fit_task.sas(kr).chi2));
            xlabel('s [Angstroem^{-1}]');
            ylabel('log(I(s))');
        end
    end
    if nr_pre > 0
        kft = 0;
        for kr = 1:length(pre_parameters)
            figure(10000+kr); clf; hold on;
            range = pre_parameters(kr).range;
            plot([1,range(2)-range(1)+1],pre_parameters(kr).max_Gamma2*[1,1],'Color',[0,0.7,0]);
            kx = 0;
            max_exp = 0;
            max_fit = 0;
            for nr_pre = range(1):range(2)
                kft = kft + 1;
                kx = kx + 1;
                pre_std = fit_task.pre(kft).exp_data(2);
                minbar = fit_task.pre(kft).exp_data(1) - pre_std;
                if minbar < 0
                    minbar = 0;
                end
                maxbar = fit_task.pre(kft).exp_data(1) + pre_std;
                if ~pre_parameters(kr).fit_rates && maxbar > 1
                    maxbar = 1;
                end
                plot([kx,kx],[minbar,maxbar],'Color',[0.5,0.5,0.5],'LineWidth',1);
                plot(kx,fit_task.pre(kft).exp_data(1),'k.','MarkerSize',14);
                plot(kx,fit_task.pre(kft).fit_data,'o','Color',[0.7,0,0],'MarkerSize',8);
                if max_exp < fit_task.pre(kft).exp_data(1)
                    max_exp = fit_task.pre(kft).exp_data(1);
                end
                if max_fit < fit_task.pre(kft).fit_data
                    max_fit = fit_task.pre(kft).fit_data;
                end
            end
            title(sprintf('PRE for labelling site %s',restraints.pre(kr).site));
            xlabel('Residue');
            set(gca,'FontSize',14);
            if pre_parameters(kr).fit_rates
                axis([1,kx,0,max([max_exp,max_fit])]);
                ylabel('\Gamma_2 [s^{-1}]');
            else
                ylabel('I_{para}/I_{dia}');
                axis([1,kx,0,1.05]);
            end
        end
        fprintf(1,'PRE %s deviation: %6.4f\n',pre_fit_type,fom_pre);
    end
end

function chi2 = sas_baseline(bsl,coeff,data)

curve = data(:,2);
errors = data(:,3);
[m,n] = size(data);
basis = data(:,4:n);
% coeff = v(1:n-3);
sim = basis*coeff' + bsl;
sc = sum(curve.*sim)/sum(sim.*sim);
chi2 = sum(((curve-sc*sim)./errors).^2)/(m-1);
% sim = basis*coeff' + v(n-3+k);
% sim = fit + bsl;
% sc = sum(curve.*sim)/sum(sim.*sim);
% chi2 = sum(((curve-sc*sim)./errors).^2)/(length(curve)-1);

function fom = sim_multi_ddr(v,fit)

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

fom = 1 - fom^(1/length(fit));

function fom = sim_multi_sas(v,fit)

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

function [fom,fit_type] = sim_multi_pre(v,predictions,parameters)

fom = 0;
n = 0;
for kr = 1:length(parameters)
    td = parameters(kr).td;
    R2dia = parameters(kr).R2dia;
    range = parameters(kr).range;
    n = n + 1 + range(2) - range(1);
    % standard deviations of input parameters 
    pre_std = predictions(range(1):range(2),2);
    % if at least one of them is zero, mean-square deviation is fitted instead of chi^2
    if sum(pre_std == 0) > 0
        pre_std = ones(size(pre_std));
        fit_type = 'mean-square';
    else
        fit_type = 'chi-square';
    end
    coeff = v/sum(v);
    Gamma2 = predictions(range(1):range(2),3:end)*coeff';
    if parameters(kr).fit_rates
        all_pre = Gamma2;
    else
        all_pre = R2dia*exp(-td*Gamma2)./(Gamma2+R2dia);
    end
    fom = fom + sum(((all_pre-predictions(range(1):range(2),1))./pre_std).^2);
end
fom = fom/n;


function record_exception(exception,logfile)

fprintf(logfile,'### ensemble fit exception %s ###\n',exception.message);