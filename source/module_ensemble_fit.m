function [entity,exceptions,failed,restraints,fit_task] = module_ensemble_fit(control,logfid)
%
% MMModel Run modelling pipeline
%
%   [exceptions,entity] = MMModel(control_file)
%   Parses control file for modelling directives and arguments and runs the
%   specified modelling pipeline
%
% INPUT
% control       control structure with fields
%               .name           'ensemblefit', unused
%               .options        options struct
%               .directives     array of struct with directives
%               .entity_output  number of the entity to be output
% logfid        file handle for log file
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
entity = [];
failed = false;

% set defaults

run_fit = true;
plot_result = false;
initial_ensemble = '';
added_conformers = '';
opt.threshold = 0.01; % conformers wih population below this threshold are discarded

opt.interactive = false; % no interactive plotting during fits
opt.blocksize = 100;
opt.skip_restraints = false;

outname = 'ensemble.ens';

restraints.ddr(1).labels{1} = '';
restraints.saxs(1).data{1} = '';
restraints.sans(1).data{1} = '';

ddr_poi = 0;
saxs_poi = 0;
sans_poi = 0;

% read restraints
for d = 1:length(control.directives)
    switch lower(control.directives(d).name)
        case 'nofit'
            run_fit = false;
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
            end
            if args < 3 % misformed restraint line
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_ensemble_fit:misformed_ddr', 'ddr restraint has less than three arguments');
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
    end
end

restraints.ddr = restraints.ddr(1:ddr_poi);
restraints.saxs = restraints.saxs(1:saxs_poi);
restraints.sans = restraints.sans(1:sans_poi);

% make file list of conformers

[initial_files,pop] = rd_ensemble_definition(initial_ensemble);
[~,~,ext] = fileparts(added_conformers);
if isempty(ext)
    added_conformers = strcat(added_conformers,'.pdb');
end
added_files = dir(added_conformers); % find all files that match the pattern

C = length(initial_files) + length(added_files);

C0 = length(initial_files); % number of conformers that must be included

if C == 0
    warnings = warnings + 1;
    exceptions{warnings} = MException('module_ensemble_fit:no_conformers', 'No conformers are specified or found.');
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
fit_task.ddr(nr).sim_distr = zeros(1,length(fit_task.r_axis));
fit_task.ddr(nr).exp_distr = [];
fit_task.ddr(nr).exp_distr_lb = [];
fit_task.ddr(nr).exp_distr_ub = [];
fit_task.ddr(nr).fit_distr = [];
for ddr_poi = 1:length(restraints.ddr)
    for kr = 1:length(restraints.ddr(ddr_poi).site1)
        if ~isempty(restraints.ddr(ddr_poi).file{kr})
            exp_data = load(restraints.ddr(ddr_poi).file{kr});
            rax_exp = 10*exp_data(:,1)'; % conversion to Angstroem
            distr_exp = exp_data(:,2)';
            lb_exp = exp_data(:,3)';
            ub_exp = exp_data(:,4)';
            fit_task.ddr(kr).exp_distr = interp1(rax_exp,distr_exp,fit_task.r_axis,'pchip',0);
            sc = 1/sum(fit_task.ddr(kr).exp_distr);
            fit_task.ddr(kr).exp_distr = sc*fit_task.ddr(kr).exp_distr;
            fit_task.ddr(kr).exp_distr_lb = sc*interp1(rax_exp,lb_exp,fit_task.r_axis,'pchip',0);
            fit_task.ddr(kr).exp_distr_ub = sc*interp1(rax_exp,ub_exp,fit_task.r_axis,'pchip',0);
        end
        if ~isempty(restraints.ddr(ddr_poi).r(kr))
            rmean = restraints.ddr(ddr_poi).r(kr);
            sigr = restraints.ddr(ddr_poi).sigr(kr);
            rax = fit_task.r_axis;
            sim_distr = exp(-(rax-rmean).^2/(2*sigr^2));
            fit_task.ddr(kr).sim_distr = sim_distr/sum(sim_distr);
        end
        fit_task.ddr(kr).distr = zeros(C,length(fit_task.r_axis));
    end
end



fit_task.ddr_valid = true(1,nr); % for logical indexing of restraints
valid_conformers = true(1,C); % for logical indexing of conformers
for c = 1:C
    % fprintf(1,'Working on conformer: %s\n',fit_task.file_list{c});
    entity = get_pdb(fit_task.file_list{c});
    block_offset = 0;
    for ddr_poi = 1:length(restraints.ddr)
        label1 = restraints.ddr(ddr_poi).labels{1};
        label2 = restraints.ddr(ddr_poi).labels{2};
        for kr = 1:length(restraints.ddr(ddr_poi).site1)
            fit_task.ddr(block_offset+kr).assignment = [ddr_poi,kr];
            site1 = restraints.ddr(ddr_poi).site1{kr};
            site2 = restraints.ddr(ddr_poi).site2{kr};
            [~,distr,entity,ddr_exceptions] = distance_distribution(entity,site1,label1,site2,label2,options);
            if ~isempty(ddr_exceptions) && ~isempty(ddr_exceptions{1})
                for ex = 1:length(ddr_exceptions)
                    warnings = warnings + 1;
                    exceptions{warnings} = ddr_exceptions{ex};
                end
            end
            if ~isempty(distr)
                fit_task.ddr(block_offset+kr).distr(c,:) = distr/sum(distr);
            else
                fprintf(logfid,'Warning: Restraint %i cannot be evaluated for conformer %i\n',kr,c);
                if opt.skip_restraints % if requested, skip restraints that cannot be evaluated for all conformers
                    fit_task.ddr_valid(block_offset+kr) = false;
                else % otherwise skip conformers for which not all restraints can be evaluated
                    valid_conformers(c) = false;
                end
            end
        end
        block_offset = block_offset + length(restraints.ddr(ddr_poi).site1);
    end
end

% restrict to the conformers for which all valid restraints could be
% evaluated
fit_task.file_list = fit_task.file_list(valid_conformers);
fit_task.pop = fit_task.pop(valid_conformers);
fit_task.pop = fit_task.pop/sum(fit_task.pop);
C = length(fit_task.pop);

fprintf(logfid,'\n');

fit_task.remaining_conformers = 1:C;
fit_task.ensemble_populations = fit_task.pop;

if run_fit
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
    
        clear fit_multi_ddr % initialize iteration counter
    
        first_run = false;
        
        if C0 > opt.blocksize - 10 % increase blocksize if necessary and raise warning 
                                   % process at least 10 additional conformers
            blocksize = C0 + round(opt.blocksize/2);
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_ensemble_fit:blocksize_increased', 'Block size was increased to %i.',blocksize);
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
        
        curr_ddr_predictions = all_ddr_predictions;
        for kr = 1:length(all_ddr_predictions) % loop over all restraints
            data = all_ddr_predictions{kr};
            data = [data(:,1:2) data(:,2+conformers)];
            curr_ddr_predictions{kr} = data;
        end

        fanonym_ddr = @(v_opt)fit_multi_ddr(v_opt,curr_ddr_predictions,opt);

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

        fprintf(logfid,'DEER restraints fit took %i h %i min %i s\n',th,tmin,ts);
        switch exitflag
            case 0
                fprintf(logfid,'Warning: Maximum number of function evaluations or iterations reached. Not converged.\n');
            case 1
                fprintf(logfid,'Convergence by mesh size.\n');
            case 2
                fprintf(logfid,'Convergence by change in population.\n');
            case 3
                fprintf(logfid,'Convergence by DEER overlap deficiency precision.\n');
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
        fprintf(logfid,'Final overlap deficiency is: %6.3f with %i conformers\n\n',fom_ddr,C0);
    end
end

% save ensemble
ofid = fopen(outname,'wt');
if ofid == -1
    warnings = warnings + 1;
    exceptions{warnings} = MException('module_ensemble_fit:ensemble_not_saved', 'Output file %s could not be written.',outname);
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

% save fit task and restraints

save([outname '.mat'],'restraints','fit_task','exceptions');

if plot_result
    dr = fit_task.r_axis(2) - fit_task.r_axis(1);
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
            %     fprintf(1,'%i) Normalization of restraint distribution: %8.4f\n',k,sum(fit_task.ddr(k).sim_distr));
            %     fprintf(1,'%i) Normalization of ensemble distribution : %8.4f\n',k,sum(fit_task.ddr(k).fit_distr));
            % fprintf(1,'%i) Overlap Gaussian                       : %8.4f\n',kr,overlap_G);
            fprintf(1,'%i) Overlap DeerNet                        : %8.4f\n',kr,overlap_E);
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
end
