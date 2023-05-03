function [entity,exceptions,failed] = module_ensemble_analysis(control,logfid)
%
% MODULE_ENSEMBLE_ANALYSIS    Analyses an ensemble
%
%   [entity,exceptions] = MODULE_ENSEMBLE_ANALYSIS(control,logfid)
%   Several functions for characterizing heterogeneous order in an ensemble
%   and for comparing two ensembles
%
% INPUT
% control       control structure with fields
%               .name           'ensembleanalysis', unused
%               .options        options struct
%               .directives     array of struct with directives
%               .entity_output  number of the entity to be output
% logfid        file handle for log file, defaults to console output
%
% OUTPUT
% entity        output entity
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
% failed        flag indicating whether the module failed
%
% SUPPORTED DIRECTIVES (an ordered list is processed)
%
% addpdb        add conformers by reading pdb files
% compare       compare two ensembles
% density       density map
% figures       format for saving figures, 'off' switches figure saving
%               off, this is the default
% flexibility   residue-specific flexibility by Ramachandran angle
%               distribution
% getens        get ensemble file
% measures      several measures of an ensemble:
%                   Rg              radius of gyration analysis
%                   width           ensemble width and density
%                   correlation     pair correlation matrix
%                   compactness     compactness analysis with Flory model
%                   chains          analysis is performed for individual
%                                   chains
%                   oriented        conformers are assumed to be already
%                                   superimposed
%                   matlab          save output as matlab files
%                   csv             save output as comma-separated value
%                                   files
% order         local order of residues with respect to the whole structure 
% sort          sort ensemble
% subsample     reduce ensemble (from a trajectory) by subsampling
% superimpose   superimpose conformers
% save          saves ensemble in a single PDB file
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% initialize empty output
exceptions = {[]};
entity = [];
warnings = 0;
failed = false;


% set defaults

save_figures = true; % default is to save figures
figure_format = 'pdf';

% default output to Matlab console if no log file identifiere was provided 
if isempty(logfid)
    logfid = 1;
end

commands = cell(1,1000); % command list
ensembles = cell(1,100); % ensemble list

cmd_poi = 0; % command pointer
ensemble_poi = 0; % entity pointer
% reorganize command line arguments
for d = 1:length(control.directives)
    clear cmd
    cmd.name = lower(control.directives(d).name);
    switch lower(control.directives(d).name)
        case {'addpdb','getens','input'}
            cmd_poi = cmd_poi + 1;
            ensemble_poi = ensemble_poi + 1;
            cmd.input = control.directives(d).options{1};
            cmd.ensemble = ensemble_poi;
            if length(control.directives(d).options) > 1 % the entity has an explicit internal name
                ensemble_descriptor.name = control.directives(d).options{2};
            else % the internal name is derived from the order of loading entities
                ensemble_descriptor.name = sprintf('E%i',ensemble_poi);
            end
            ensembles{ensemble_poi} = ensemble_descriptor;
            commands{cmd_poi} = cmd;
        case 'figures'
            cmd_poi = cmd_poi + 1;
            cmd.extension = control.directives(d).options{1};
            commands{cmd_poi} = cmd;
        case 'flexibility'
            cmd_poi = cmd_poi + 1;
            cmd.fname = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % a selected entity is analyzed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; % flexibility analysis is performed for current entity
            end
            commands{cmd_poi} = cmd;
        case 'order'
            cmd_poi = cmd_poi + 1;
            cmd.fname = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % a selected entity is analyzed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; % Local order analysis is performed for current entity
            end
            commands{cmd_poi} = cmd;
        case 'compare'
            cmd_poi = cmd_poi + 1;
            cmd.entity1 = control.directives(d).options{1};
            cmd.entity2 = control.directives(d).options{2};
            cmd.address = '';
            cmd.resolved = false;
            if length(control.directives(d).options) > 2 % chain and possibly range given
                cmd.address = control.directives(d).options{3};
            end
            if length(control.directives(d).options) > 3 
                if strcmpi(control.directives(d).options{4},'resolved')
                    cmd.resolved = true;
                end
            end
            commands{cmd_poi} = cmd;
        case 'subsample'
            cmd_poi = cmd_poi + 1;
            cmd.stride = str2double(control.directives(d).options{1});
            if length(control.directives(d).options) > 1 % the entity has an explicit internal name
                cmd.entity = control.directives(d).options{2};
            else % use current entity
                cmd.entity = '.';
            end
            if length(control.directives(d).options) > 2 % the entity has an explicit internal name
                cmd.output = control.directives(d).options{3};
            else % use current entity
                cmd.output = '.';
            end
            commands{cmd_poi} = cmd;
        case 'sort'
            cmd_poi = cmd_poi + 1;
            cmd.outname = control.directives(d).options{1};
            cmd.oriented = false;
            cmd.drms = false;
            if length(control.directives(d).options) > 1 % a selected entity is analyzed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; 
            end
            if length(control.directives(d).options) > 2 % further directive
                if strcmpi(control.directives(d).options{3},'drms')
                    cmd.drms = true;
                end
                if strcmpi(control.directives(d).options{3},'oriented')
                    cmd.oriented = true;
                end
            end
            commands{cmd_poi} = cmd;
        case 'save'
            cmd_poi = cmd_poi + 1;
            cmd.outname = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % a selected entity is analyzed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; 
            end
            commands{cmd_poi} = cmd;
        case 'measures'
            cmd_poi = cmd_poi + 1;
            cmd.outname = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % a selected entity is analyzed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; 
            end
            cmd.options.chain_mode = false;
            cmd.options.Rg = false;
            cmd.options.pair_rmsd = false;
            cmd.options.pair_drms = false;
            cmd.options.superimpose = true;
            cmd.options.sorted = false;
            cmd.options.pair_corr = false;
            cmd.options.compactness = false;
            cmd.options.matlab = false;
            cmd.options.csv = false;
            cmd.options.chain = '';
            cmd.options.range = [];
            if length(control.directives(d).options) > 2 % chain and possibly range given
                cmd.options.chain_mode = true;
                [chain,range] = split_chain_range(control.directives(d).options{3});
                cmd.options.chain = chain;
                cmd.options.range = range;
            end
            [n,~] = size(control.directives(d).block);
            % set requested options
            for k = 1:n
                switch lower(control.directives(d).block{k,1})
                    case 'rg'
                        cmd.options.Rg = true;
                    case 'width'
                        cmd.options.pair_rmsd = true;
                    case 'drms'
                        cmd.options.pair_drms = true;
                    case 'sort'
                        cmd.options.sorted = true;
                    case 'correlation'
                        cmd.options.pair_corr = true;
                    case 'compactness'
                        cmd.options.compactness = true;
                    case 'oriented'
                        cmd.options.superimpose = false;
                    case 'matlab'
                        cmd.options.matlab = true;
                    case 'csv'
                        cmd.options.csv = true;
                end
            end
            commands{cmd_poi} = cmd;            
        case 'density'
            cmd_poi = cmd_poi + 1;
            cmd.outname = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % a selected entity is analyzed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; 
            end
            cmd.address = '(*)'; % by default, all chains are selected
            if length(control.directives(d).options) > 2 % chain and possibly range given
                cmd.address = control.directives(d).options{3};
            end
            cmd.resolution = 1;
            if length(control.directives(d).options) > 3 
                cmd.resolution = str2double(control.directives(d).options{4});
            end
            commands{cmd_poi} = cmd;            
        case 'property'
            cmd_poi = cmd_poi + 1;
            cmd.outname = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % a selected entity is analyzed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; 
            end
            cmd.address = '(*)'; % by default, all chains are selected
            if length(control.directives(d).options) > 2 % chain and possibly range given
                cmd.address = control.directives(d).options{3};
            end
            cmd.resolution = 1;
            if length(control.directives(d).options) > 3 
                cmd.resolution = str2double(control.directives(d).options{4});
            end
            cmd.property = 'electrostatic';
            if length(control.directives(d).options) > 4 
                cmd.property = control.directives(d).options{5};
            end
            cmd.pH = 7;
            if length(control.directives(d).options) > 5 
                cmd.pH = str2double(control.directives(d).options{6});
            end
            cmd.I = 0.150;
            if length(control.directives(d).options) > 7 
                cmd.I = str2double(control.directives(d).options{7});
            end
            commands{cmd_poi} = cmd;            
        case 'coulomb'
            cmd_poi = cmd_poi + 1;
            cmd.outname = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % a selected entity is analyzed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; 
            end
            cmd.address = '(A)'; % by default, chain A is selected
            if length(control.directives(d).options) > 2 % chain and possibly range given
                cmd.address = control.directives(d).options{3};
            end
            cmd.aa1 = 'Arg';
            if length(control.directives(d).options) > 3 
                cmd.aa1 = control.directives(d).options{4};
            end
            cmd.aa2 = 'Glu';
            if length(control.directives(d).options) > 4 
                cmd.aa2 = control.directives(d).options{5};
            end
            cmd.pH = 7;
            if length(control.directives(d).options) > 5 
                cmd.pH = str2double(control.directives(d).options{6});
            end
            cmd.I = 0.150;
            if length(control.directives(d).options) > 6 
                cmd.I = str2double(control.directives(d).options{7});
            end
            cmd.maxscale = [];
            if length(control.directives(d).options) > 7 
                cmd.maxscale = str2double(control.directives(d).options{8});
            end
            commands{cmd_poi} = cmd;            
        case 'superimpose'
            cmd_poi = cmd_poi + 1;
            cmd.outname = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % a selected entity is analyzed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; % superposition is performed for current entity
            end
            cmd.address = '';
            if length(control.directives(d).options) > 2 % chain and possibly range given
                cmd.address = control.directives(d).options{3};
            end
            cmd.template = '';
            if length(control.directives(d).options) > 3 % superposition to a template
                cmd.template = control.directives(d).options{4};
            end
            cmd.template_address = '';
            if length(control.directives(d).options) > 4 % chain and possibly range given for template
                cmd.template_address = control.directives(d).options{5};
            end
            cmd.central = false;
            if length(control.directives(d).options) > 5 && strcmpi(control.directives(d).options{6},'central')
                cmd.central = true;
            end
            commands{cmd_poi} = cmd;    
        case 'inertiaframe'
            cmd_poi = cmd_poi + 1;
            cmd.outname = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % a selected entity is analyzed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; % superposition is performed for current entity
            end
            cmd.address = '';
            if length(control.directives(d).options) > 2 % chain and possibly range given for template
                cmd.address = control.directives(d).options{3};
            end
            commands{cmd_poi} = cmd;    
        case 'asphericity'
            cmd_poi = cmd_poi + 1;
            if ~isempty(control.directives(d).options) % a selected entity is analyzed
                cmd.entity = control.directives(d).options{1};
            else
                cmd.entity = '.'; % asphericity is computed for current entity
            end
            cmd.address = '';
            if length(control.directives(d).options) > 1 % chain and possibly range given
                cmd.address = control.directives(d).options{1};
            end
            commands{cmd_poi} = cmd;    
        otherwise
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_ensembleanalysis:unknown_directive',...
                'directive %s is unknown',lower(control.directives(d).name));
            record_exception(exceptions{warnings},logfid);
    end
end

ensembles = ensembles(1:ensemble_poi);
commands = commands(1:cmd_poi);

fprintf(logfid,'\n%i commands will be executed on %i ensembles\n\n',cmd_poi,ensemble_poi);

% run the command list

for c = 1:cmd_poi
    cmd = commands{c};
    switch cmd.name
        case 'addpdb'
            ensemble_poi = cmd.ensemble;
            ensemble_descriptor = ensembles{ensemble_poi};
            ensemble_name = ensemble_descriptor.name;
            [~,~,ext] = fileparts(cmd.input);
            if isempty(ext)
                cmd.input = strcat(cmd.input,'.pdb');
            end
            added_files = dir(cmd.input); % find all files that match the pattern
            [entity,exceptions] = get_pdb(added_files(1).name);
            if ~isempty(exceptions) && ~isempty(exceptions{1})
                warnings = warnings +1;
                exceptions{warnings} = MException('module_ensembleanalysis:file_does_not_exist',...
                    'PDB file %s could not be opened',added_files(1).name);
                record_exception(exceptions{warnings},logfid);
                return
            end
            if length(added_files) > 1
                for cf = 2:length(added_files)
                    [entity,exceptions] = get_pdb(added_files(cf).name,[],entity);
                    if ~isempty(exceptions) && ~isempty(exceptions{1})
                        warnings = warnings +1;
                        exceptions{warnings} = MException('module_ensembleanalysis:file_does_not_exist',...
                            'PDB file %s could not be opened',added_files(cf).name);
                        record_exception(exceptions{warnings},logfid);
                        return
                    end
                end
                entity.populations = ones(1,length(added_files))/length(added_files);
            end
            ensemble_descriptor.entity = entity;
            fprintf(logfid,'\nCurrent ensemble is: %s\n',ensemble_name);
            ensembles = store_ensemble(ensemble_name,entity,ensembles);
        case {'getens','input'}
            ensemble_poi = cmd.ensemble;
            ensemble_descriptor = ensembles{ensemble_poi};
            ensemble_name = ensemble_descriptor.name;
            entity = get_ensemble(cmd.input);
            fprintf(logfid,'\nCurrent ensemble is: %s\n',ensemble_name);
            ensembles = store_ensemble(ensemble_name,entity,ensembles);
        case 'figures'
            if strcmpi(cmd.extension,'off')
                save_figures = false;
            else
                save_figures = true;
                figure_format = cmd.extension;
            end
        case 'flexibility'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'flexibility cannot be analysed for entity %s, since entity is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            if save_figures
                flexibility_analysis(c_entity,cmd.fname,figure_format);
            else
                flexibility_analysis(c_entity,cmd.fname);
            end
        case 'order'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'Local order analysis cannot be performed for entity %s, since entity is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            if save_figures
                local_order(c_entity,cmd.fname,figure_format);
            else
                local_order(c_entity,cmd.fname);
            end
        case 'compare'
           entity1 = retrieve_ensemble(cmd.entity1,ensembles,logfid);
           if isempty(entity1)
               warnings = warnings +1;
               exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                   'tried to comparison with entity %s, which is unknown',cmd.entity1);
               record_exception(exceptions{warnings},logfid);
               return
           end
           if ~strcmpi(cmd.entity2,'-self')
               entity2 = retrieve_ensemble(cmd.entity2,ensembles,logfid);
               if isempty(entity2)
                   warnings = warnings +1;
                   exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                       'tried to comparison with entity %s, which is unknown',cmd.entity2);
                   record_exception(exceptions{warnings},logfid);
                   return
               end
           else
               entity2 = [];
           end
           if cmd.resolved
               [chain,range] = split_chain_range(cmd.address);
               residues = range(1):range(2);
               overlaps = zeros(size(residues));
               for r = residues
                   address = sprintf('(%s)%i',chain,r);
                   fprintf(1,'Computing overlap for %s\n',address);
                   overlaps(r-residues(1)+1) = density_overlap(entity1,entity2,address);
               end
               h = figure;
               plot(residues,overlaps,'.','MarkerSize',14,'Color',[0.25,0.25,0.25]);
               xlabel('Residue number');
               ylabel('Overlap');
               axis([range(1),range(2),0,1]);
               title(sprintf('Overlap between ensembles %s and %s',cmd.entity1,cmd.entity2));
               if save_figures
                   figname = sprintf('overlap_%s_%s.%s',cmd.entity1,cmd.entity2,figure_format);
                   saveas(h,figname);
               end
               datname = sprintf('overlap_%s_%s.mat',cmd.entity1,cmd.entity2);
               save(datname,'residues','overlaps');
               data = [residues' overlaps'];
               datname = sprintf('overlap_%s_%s.csv',cmd.entity1,cmd.entity2);
               writematrix(data,datname);
           else
               overlap = density_overlap(entity1,entity2,cmd.address);
               fprintf(logfid,'Ensemble density overlap between %s and %s is %6.3f\n',cmd.entity1,cmd.entity2,overlap);
           end
        case 'subsample'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'tried to subsample entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            C = length(c_entity.populations);
            selection = '{1';
            conf = 1;
            while conf + cmd.stride <= C
                conf = conf + cmd.stride;
                selection = sprintf('%s,%i',selection,conf);
            end
            selection = sprintf('%s}(*)',selection);
            c_entity = select(c_entity,selection,true);
            clear save_options
            save_options.selected = true;
            save_options.pop = true;
            put_pdb(c_entity,'temporary.pdb',save_options);
            c_entity = get_pdb('temporary.pdb');
            delete('temporary.pdb');
            if strcmp(cmd.output,'.')
                entity = c_entity;
            end
            ensembles = store_ensemble(cmd.output,c_entity,ensembles);
        case 'sort'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'tried to subsample entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            if cmd.drms
                [pair_rmsd,pop,exceptions0] = pair_drms_matrix(c_entity);
            elseif cmd.oriented
                [pair_rmsd,pop,exceptions0] = pair_rmsd_matrix_oriented(c_entity);
            else
                [pair_rmsd,pop,exceptions0] = pair_rmsd_matrix(c_entity);
            end
            if ~isempty(exceptions0{1})
                for k = 1:exceptions0
                    warnings = warnings + 1;
                    exceptions{warnings} = exceptions0{k};
                end
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_ensembleanalysis:backbone_retrieval_failed',...
                    'sorting of ensemble %s failed since backbone could not be retrieved',cmd.entity);
                record_exception(exceptions{warnings},logfid);
                return
            end                
            [pair_rmsd,ordering,~,cluster_sizes,cluster_pop] = cluster_sorting(pair_rmsd,pop);
            [pname,fname,~] = fileparts(cmd.outname);
            basname = fullfile(pname,fname);
            ensemble_name = strcat(basname,'.ens');
            ens_fid = fopen(ensemble_name,'wt');
            fprintf(ens_fid,'%% Sorted ensemble %s by MMMx\n',basname);
            poi = 0;
            for clust = 1:length(cluster_pop)
                fprintf(logfid,'\nCluster %i with population %6.4f comprises the following %i conformers:\n',clust,cluster_pop(clust),cluster_sizes(clust));
                for conf = 1:cluster_sizes(clust)
                    poi = poi + 1;
                    oname = sprintf('%s_m%i.pdb',basname,poi);
                    fprintf(logfid,'%s (conformer %i in the input ensemble)\n',oname,ordering(poi));
                    clear save_options
                    save_options.order = ordering(poi);
                    exceptions = put_pdb(c_entity,oname,save_options);
                    fprintf(ens_fid,'%s  %8.6f %% cluster %i\n',oname,pop(ordering(poi)),clust);
                end
            end
            fclose(ens_fid);
            if cmd.drms
                h = plot_pair_rmsd(pair_rmsd,false,true);
                if save_figures
                    figname = sprintf('pair_drms_sorting_%s.%s',basname,figure_format);
                    saveas(h,figname);
                end
            else
                h = plot_pair_rmsd(pair_rmsd);
                if save_figures
                    figname = sprintf('pair_rmsd_sorting_%s.%s',basname,figure_format);
                    saveas(h,figname);
                end
            end
            c_entity = get_ensemble(ensemble_name);
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            ensembles = store_ensemble(cmd.entity,c_entity,ensembles);
        case 'save'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'tried to subsample entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            save_options.pop = true;
            exceptions0 = put_pdb(c_entity,cmd.outname,save_options);
            for exci = 1:length(exceptions0)
                if ~isempty(exceptions0{exci})
                    warnings = warnings +1;
                    exceptions{warnings} = exceptions0{exci};
                    record_exception(exceptions{warnings},logfid);
                end
            end
        case 'measures'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'tried to subsample entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            [pname,fname,~] = fileparts(cmd.outname);
            basname = fullfile(pname,fname);
            if cmd.options.chain_mode
                [backbones,pop,exceptions0] = get_backbones_ensemble(c_entity,cmd.options.chain,cmd.options.range);
            else
                [backbones,pop,exceptions0] = get_backbones_ensemble(c_entity);
            end
            if ~isempty(exceptions0{1})
                for k = 1:exceptions0
                    warnings = warnings + 1;
                    exceptions{warnings} = exceptions0{k};
                end
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_ensembleanalysis:backbone_retrieval_failed',...
                    'measures for ensemble %s cannot be computed since backbone could not be retrieved',cmd.entity);
                record_exception(exceptions{warnings},logfid);
                return
            end
            [measures,correlations] = analyze_ensemble(backbones,pop,cmd.options);
            parts = fieldnames(measures);
            for p = 1:length(parts)
                part = parts{p};
                if strcmp(part,'all')
                    name = 'whole structure';
                else
                    name = sprintf('chain %s',part);
                end
                if cmd.options.Rg
                    fprintf(logfid,'\nRadius of gyration analysis for %s (%s)\n',name,part);
                    fprintf(logfid,'   Rg = %4.1f %s with standard deviation %4.1f %s\n',...
                        measures.(part).Rg,char(197),measures.(part).Rg_std,char(197));
                end
                if cmd.options.pair_rmsd
                    fprintf(logfid,'\nEnsemble width and density for %s (%s)\n',name,part);
                    fprintf(logfid,'   width = %4.1f %s; density %4.1f %s\n',...
                        measures.(part).width,char(197),measures.(part).density,char(197));
                    h = plot_pair_rmsd(measures.(part).pair_rmsd,cmd.options.superimpose);
                    sum_msq = sum(measures.(part).pair_rmsd.^2);
                    [C,~] = size(measures.(part).pair_rmsd);
                    [msq,conf] = min(sum_msq);
                    rmsd = sqrt(msq/C);
                    fprintf(logfid,'\nCentral conformer %i has rmsd distance of %4.1f %s from other conformers\n',...
                        conf,rmsd,char(197));
                    if save_figures
                        figname = sprintf('pair_rmsd_%s_%s.%s',basname,part,figure_format);
                        saveas(h,figname);
                    end
                    if cmd.options.matlab
                        pair_rmsd = measures.(part).pair_rmsd;
                        datname = sprintf('pair_rmsd_%s_%s.mat',basname,part);
                        save(datname,'pair_rmsd');
                    end
                    if cmd.options.csv
                        pair_rmsd = measures.(part).pair_rmsd;
                        datname = sprintf('pair_rmsd_%s_%s.csv',basname,part);
                        writematrix(pair_rmsd,datname);
                    end
                end
                if cmd.options.pair_drms
                    h = plot_pair_rmsd(measures.(part).pair_drms,cmd.options.superimpose,cmd.options.pair_drms);
                    sum_msq = sum(measures.(part).pair_drms.^2);
                    [C,~] = size(measures.(part).pair_drms);
                    [msq,conf] = min(sum_msq);
                    drms = sqrt(msq/C);
                    fprintf(logfid,'\nCentral conformer %i has drms of %4.1f %s from other conformers\n',...
                        conf,drms,char(197));
                    if save_figures
                        figname = sprintf('pair_drms_%s_%s.%s',basname,part,figure_format);
                        saveas(h,figname);
                    end
                    if cmd.options.matlab
                        pair_drms = measures.(part).pair_drms;
                        datname = sprintf('pair_drms_%s_%s.mat',basname,part);
                        save(datname,'pair_drms');
                    end
                    if cmd.options.csv
                        pair_drms = measures.(part).pair_drms;
                        datname = sprintf('pair_drms_%s_%s.csv',basname,part);
                        writematrix(pair_drms,datname);
                    end
                end
                if cmd.options.pair_corr
                    h = plot_pair_corr(correlations.(part).pair_axis,correlations.(part).pair_corr);
                    if save_figures
                        figname = sprintf('residue_pair_correlation_%s_%s.%s',basname,part,figure_format);
                        saveas(h,figname);
                    end
                    h = plot_pair_corr(correlations.(part).pair_axis,correlations.(part).pair_corr,correlations.(part).dmat);
                    if save_figures
                        figname = sprintf('residue_pair_correlation_normalized_%s_%s.%s',basname,part,figure_format);
                        saveas(h,figname);
                    end
                    if cmd.options.matlab
                        pair_corr = correlations.(part).pair_corr;
                        dmat = correlations.(part).dmat;
                        datname = sprintf('residue_pair_correlation_%s_%s.mat',basname,part);
                        save(datname,'pair_corr','dmat');
                    end
                    if cmd.options.csv
                        pair_corr = correlations.(part).pair_corr;
                        dmat = correlations.(part).dmat;
                        datname = sprintf('residue_pair_correlation_%s_%s.csv',basname,part);
                        writematrix(pair_corr,datname);
                        datname = sprintf('distance_matrix_%s_%s.csv',basname,part);
                        writematrix(dmat,datname);
                    end
                end
                if cmd.options.compactness
                    fprintf(logfid,'\nCompactness analysis for %s (%s)\n',name,part);
                    offset = 0;
                    if ~isempty(cmd.options.range)
                        offset = cmd.options.range(1); 
                    end                    
                    compactness = correlations.(part).compact;
                    h = plot_compactness(compactness,'Compactness matrix',[-1,1],offset);
                    if save_figures
                        figname = sprintf('compactness_matrix_%s_%s.%s',basname,part,figure_format);
                        saveas(h,figname);
                    end
                    proximity = correlations.(part).proximity;
                    h = plot_compactness(proximity,'Proximity matrix',[-1,1],offset);                    
                    if save_figures
                        figname = sprintf('proximity_matrix_%s_%s.%s',basname,part,figure_format);
                        saveas(h,figname);
                    end
                    Rdev = correlations.(part).Rdev;
                    extent = max([abs(min(min(Rdev))),max(max(Rdev))]);
                    h = plot_compactness(Rdev,'Segment length deviation',1.05*[-extent,extent],offset);
                    if save_figures
                        figname = sprintf('segment_length_deviation_%s_%s.%s',basname,part,figure_format);
                        saveas(h,figname);
                    end
                    h = plot_segments(measures.(part));
                    if save_figures
                        figname = sprintf('segment_length_distribution_%s_%s.%s',basname,part,figure_format);
                        saveas(h,figname);
                    end
                    seg_lengths = measures.(part).seg_lengths;
                    R2fct = measures.(part).R0_seglen*seg_lengths.^measures.(part).nu_seglen;
                    mean_R2 = measures.(part).mean_R2;                    
                    min_R2 = measures.(part).min_R2;                    
                    max_R2 = measures.(part).max_R2;
                    fprintf(logfid,'   Random coil fit\n   R0 = %5.2f; nu %5.3f\n',...
                        measures.(part).R0_seglen,measures.(part).nu_seglen);
                    segment_data = [seg_lengths' mean_R2' min_R2' max_R2' R2fct']; 
                    if cmd.options.matlab
                        datname = sprintf('compactness_analysis_%s_%s.mat',basname,part);
                        save(datname,'compactness','proximity','Rdev','seg_lengths','mean_R2','min_R2','max_R2','R2fct');
                    end
                    if cmd.options.csv
                        datname = sprintf('compactness_matrix_%s_%s.csv',basname,part);
                        writematrix(compactness,datname);
                        datname = sprintf('proximity_matrix_%s_%s.csv',basname,part);
                        writematrix(proximity,datname);
                        datname = sprintf('segment_length_deviation_%s_%s.csv',basname,part);
                        writematrix(Rdev,datname);
                        datname = sprintf('segment_length_distribution_%s_%s.csv',basname,part);
                        writematrix(segment_data,datname);
                    end
                end
            end
        case 'density'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'tried to subsample entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            make_density(c_entity,cmd.outname,cmd.address,cmd.resolution);
        case 'property'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'tried to compute property cube for entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            property_cube(c_entity,cmd.outname,cmd.property,cmd.address,cmd.resolution,cmd.pH,cmd.I)
        case 'coulomb'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'tried to compute Coulomb pairs for entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            [pname,fname,~] = fileparts(cmd.outname);
            basname = fullfile(pname,fname);
            [pairs,coulomb,~,resnums] = coulomb_pairs(c_entity,cmd.address,cmd.aa1,cmd.aa2,cmd.pH,cmd.I);
            coulomb_matrix = [pairs coulomb];
            h = plot_interaction(pairs,coulomb,resnums,cmd.maxscale);
            if save_figures
                figname = sprintf('coulomb_interaction_%s_%s_%s.%s',basname,cmd.aa1,cmd.aa2,figure_format);
                exportgraphics(h,figname);
            end
            datname = sprintf('coulomb_matrix_%s_%s_%s.csv',cmd.aa1,cmd.aa2,basname);
            writematrix(coulomb_matrix,datname);
        case 'superimpose'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'tried to superimpose entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            if isempty(cmd.template)
                t_entity = [];
            else
                t_entity = retrieve_ensemble(cmd.template,ensembles,logfid);
                if isempty(t_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'tried to superimpose onto template %s, which is unknown',cmd.template);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            [pname,fname,~] = fileparts(cmd.outname);
            basname = fullfile(pname,fname);
            selected = '';
            if ~isempty(cmd.address)
                if isempty(cmd.template_address)
                    selected = cmd.address;
                else
                    clear selected
                    selected{1} = cmd.address;
                    selected{2} = cmd.template_address;
                end
            end
            [c_entity,~,exceptions0] = superimpose_ensemble(c_entity,selected,t_entity,cmd.central);
            for exci = 1:length(exceptions0)
                if ~isempty(exceptions0{exci})
                    warnings = warnings +1;
                    exceptions{warnings} = exceptions0{exci};
                    record_exception(exceptions{warnings},logfid);
                end
            end
            pop = c_entity.populations;
            ensemble_name = strcat(basname,'.ens');
            ens_fid = fopen(ensemble_name,'wt');
            for conf = 1:length(pop)
                oname = sprintf('%s_m%i.pdb',basname,conf);
                clear save_options
                save_options.order = conf;
                exceptions = put_pdb(c_entity,oname,save_options);
                fprintf(ens_fid,'%s  %8.6f\n',oname,pop(conf));
            end
            fclose(ens_fid);
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            ensembles = store_ensemble(cmd.entity,c_entity,ensembles);
        case 'inertiaframe'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'tried to transform entity %s, which is unknown, to inertia frame',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            [pname,fname,~] = fileparts(cmd.outname);
            basname = fullfile(pname,fname);
            selected = '';
            if ~isempty(cmd.address)
                selected = cmd.address;
            end
            c_entity = inertia_frame(c_entity,selected);
            pop = c_entity.populations;
            ensemble_name = strcat(basname,'.ens');
            ens_fid = fopen(ensemble_name,'wt');
            for conf = 1:length(pop)
                oname = sprintf('%s_m%i.pdb',basname,conf);
                clear save_options
                save_options.order = conf;
                exceptions = put_pdb(c_entity,oname,save_options);
                fprintf(ens_fid,'%s  %8.6f\n',oname,pop(conf));
            end
            fclose(ens_fid);
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            ensembles = store_ensemble(cmd.entity,c_entity,ensembles);
        case 'asphericity'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'tried to compute asphericity of entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            selected = '';
            if ~isempty(cmd.address)
                selected = cmd.address;
            end
            c_entity = asphericity(c_entity,selected);
            pop = c_entity.populations;
            h = plot_correlation(c_entity.Rg_c,c_entity.asphericity_c,pop);
            ha = h.CurrentAxes;
            ha.Title.String = sprintf('Mean asphericity %5.3f',c_entity.asphericity);
            ha.XLabel.String = sprintf('R_g (%c)',char(197));
            ha.YLabel.String = 'Asphericity';
            if save_figures
                figname = sprintf('Asphericity_%s.%s',cmd.entity,figure_format);
                saveas(h,figname);
            end
            asphericity_matrix = [c_entity.populations c_entity.Rg_c c_entity.asphericity_c];
            datname = sprintf('asphericity_%s.csv',cmd.entity);
            writematrix(asphericity_matrix,datname);
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            ensembles = store_ensemble(cmd.entity,c_entity,ensembles);
    end
end

function entity = retrieve_ensemble(name,ensembles,logfid)

entity = [];
e = 0;
while isempty(entity) && e < length(ensembles)
    e = e + 1;
    ensemble_descriptor = ensembles{e};
    if strcmp(name,ensemble_descriptor.name)
        entity = ensemble_descriptor.entity;
    end
end
if isempty(entity)
    fprintf(logfid,'Ensemble %s has not been loaded.\n',name);
end

function ensembles = store_ensemble(name,entity,ensembles)

found = false;
e = 0;
while ~found && e < length(ensembles)
    e = e + 1;
    ensemble_descriptor = ensembles{e};
    if strcmp(name,ensemble_descriptor.name)
        found = true;
        ensemble_descriptor.entity = entity;
        ensembles{e} = ensemble_descriptor;
    end
end

function record_exception(exception,logfid)

fprintf(logfid,'### ensembleanalysis exception: %s ###\n',exception.message);

function [chain,range] = split_chain_range(address)

range = [];

poia = strfind(address,'(');
poie = strfind(address,')');
if isempty(poie)
    chain = '';
else
    chain = address(poia+1:poie-1);
    address = address(poie+1:end);
end
if strcmp(chain,'*')
    range = [1,1e6];
    return
end
residues = split(address,'-');
if ~isempty(residues{1})
    range(1) = str2double(residues{1});
end
if ~isempty(residues{2})
    range(2) = str2double(residues{2});
end
if length(range) == 1
    range(2) = range(1);
end


function h = plot_pair_rmsd(pair_rmsd,superimposed,drms)

if ~exist('superimposed','var') || isempty(superimposed)
    superimposed = true;
end

if ~exist('drms','var') || isempty(drms)
    drms = false;
end

if drms
    superimposed = false;
end

h = figure;

[n1,n2] = size(pair_rmsd);
plot(1,1,'k.');
plot(n1,n2,'k.');
image(pair_rmsd,'CDataMapping','scaled','ButtonDownFcn',@(hObject,eventdata,handles)axes_rmsd_ButtonDownFcn);
curr_axis = gca;
curr_axis.YDir = 'normal';
colorbar;
axis tight
xlabel('Conformer number');
ylabel('Conformer number');
if drms
    title('Pair distance root mean square deviation');
elseif superimposed
    title('Pair rmsd upon optimal superposition');
else
    title('Pair rmsd for original orientation');
end
axis equal

function h = plot_pair_corr(pair_axis,pair_corr,dmat)

h = figure;
[n1,n2] = size(pair_corr);
plot(1,1,'k.');
plot(n1,n2,'k.');

if ~exist('dmat','var') || isempty(dmat)
    image(pair_corr,'CDataMapping','scaled');
    title('Pair correlation');
else
    image(pair_corr./dmat,'CDataMapping','scaled');
    title('Pair correlation (normalized)');
end

curr_axis = gca;
curr_axis.YDir = 'normal';
step =floor(n1/5);
curr_axis.XTick = 1:step:n2;
curr_axis.YTick = 1:step:n1;
XTickLabel = cell(1,10);
poi = 0;
for k = 1:step:n2
    poi = poi + 1;
    XTickLabel{poi} = pair_axis{k};
end
XTickLabel = XTickLabel(1:poi);
curr_axis.XTickLabel = XTickLabel;

YTickLabel = cell(1,10);
poi = 0;
for k = 1:step:n2
    poi = poi + 1;
    YTickLabel{poi} = pair_axis{k};
end
YTickLabel = YTickLabel(1:poi);
curr_axis.YTickLabel = YTickLabel;
colorbar;
axis tight
xlabel('Residue number');
ylabel('Residue number');
axis equal


function h = plot_compactness(compactness,type,range,offset)

if ~exist('type','var') || isempty(type)
    type = 'Compactness matrix';
end

if ~exist('offset','var') || isempty(offset)
    offset = 0;
end

% make blue-red color map
mymap = ones(51,3);
for k = 1:25
    mymap(k,2:3) = k/26*[1,1];
    mymap(k+26,1:2) = (25-k)/25*[1,1];
end
mymap = flipud(mymap);

h = figure;
image(offset,offset,compactness,'CDataMapping','scaled');
hold on;
curr_axis = gca;
set(curr_axis,'CLim',range);
set(curr_axis,'YDir','normal');
colorbar;
axis tight
xlabel('Residue number');
ylabel('Residue number');
axis equal
colormap(mymap)
title(type);

function h = plot_segments(segments)

h = figure;

kaxis = 1:max(segments.seg_len);

% h1 = plot(segments.seg_len,segments.all_R2,'k.','MarkerSize',8);
h1 = fill([segments.seg_lengths, fliplr(segments.seg_lengths)],[segments.max_R2, fliplr(segments.min_R2)],0.75*[1,1,1],'LineStyle','none');
hold on
R2fct = segments.R0_seglen*kaxis.^segments.nu_seglen;
h2 = plot(kaxis,segments.mean_R2,'Color',[0,0.75,0],'LineWidth',2.5);
h3 = plot(kaxis,R2fct,'Color',[0.8,0,0],'LineWidth',2.5);
xlabel('Segment sequence length k');
ylabel(sprintf('<R^{2}>^{1/2} [%s]',char(197)));
legend([h1,h2,h3],'segment length distribution','mean value',sprintf('random coil %5.3f k^{%5.3f}',segments.R0_seglen,segments.nu_seglen),'Location','southeast');
axis([min(kaxis)-1,max(kaxis)+1,0,1.05*max(segments.max_R2)]);

function h = plot_correlation(x,y,pop)

h = figure;
hold on
for k = 1:length(pop)
    MS = 24*round(sqrt(pop(k)/max(pop)));
    if MS > 0
        plot(x(k),y(k),'k.','MarkerSize',MS);
    end
end

function h = plot_interaction(pairs,interaction,resnums,maxscale)

interaction = abs(interaction);
max_interaction = max(interaction);
if isempty(maxscale)
    maxscale = max_interaction;
end
h = figure;
hold on;
for k = 1:length(interaction)
    col = color_grade(interaction(k),maxscale);
    plot(pairs(k,1),pairs(k,2),'.','MarkerSize',12,'Color',col);
    plot(pairs(k,2),pairs(k,1),'.','MarkerSize',12,'Color',col);
end
curr_axis = gca;
set(curr_axis,'Color','k');
xlabel('Residue number');
ylabel('Residue number');
axis equal
axis([min(resnums),max(resnums),min(resnums),max(resnums)]);
title(sprintf('Mean Coulomb interaction with maximum %4.1f K',maxscale));



function color = color_grade(p,P)
% color grade in a hot map from 0 to P corresponding to value p  

if p < 0
    p = 0;
end
if p > P
    p = P;
end
k = 1 + round(100*p/P);

colmap = hot(101);
color = colmap(k,:);