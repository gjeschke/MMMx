function [entity,exceptions,failed] = module_experiment_design(control,logfid)
%
% MODULE_EXPERIMENT_DESIGN    Various functions for designing spin-label
%                             based experiments
%
%   [entity,exceptions,failed] = MODULE_EXPERIMENT_DESIGN(control,logfid)
%   Computes spin label rotamer distributions, distance distributions, and
%   determines suitable site pairs for modelling projects
%
% INPUT
% control       control structure with fields
%               .name           'experimentdesign', unused
%               .options        options struct
%               .directives     array of struct with directives
%               .entity_output  number of the entity to be output
% logfid        file handle for log file, defaults to console output
%
% OUTPUT
% entity        output entity
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
% failed        flag indicating whether the module failed, failure is only
%               signalled if the last entity could not be loaded
%
% SUPPORTED DIRECTIVES (an ordered list is processed)
%
% addpdb        add conformers by reading pdb files
% getens        get ensemble file
% expand        expand rigid-body arrangement 
% import        import entity directly from PDB server
% sitescan      spin-labelling site scan
% pairlist      list of generally suitable site pairs
% plot          plot distrubtions (with basis file name for saving plots)
% distributions compute and save distance distributions
% enmPairs      predicts favorable site pairs for elastic network modeling
%               of a conformation transition
% rbReference   selects reference points in rigid bodies for a project that
%               aims to use Rigi
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% initialize empty output
exceptions = {[]};
entity = [];
warnings = 0;
failed = false;

% set defaults

save_figures = false; % default is to not save figures
figure_format = 'pdf';

min_rotamers = 1; % minimum number of rotamers fo suitable labelling site
min_partition_function = 0.1; % minimum partition function
min_distance = 20; % minimum mean distance for suitable site pair
max_distance = 60; % maximum mean distance for suitable site pair
residue_types = 'CILMSTV';

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
        case {'addpdb','getens','expand','import','input','getalphafold'}
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
        case 'plot'
            cmd_poi = cmd_poi + 1;
            if ~isempty(control.directives(d).options)
                cmd.fname = control.directives(d).options{1};
            else
                cmd.fname = '';
            end
            if length(control.directives(d).options) > 1
                cmd.extension = control.directives(d).options{2};
            else
                cmd.extension = figure_format;
            end
            commands{cmd_poi} = cmd;
        case {'sitescan','sitescan!'}
            cmd_poi = cmd_poi + 1;
            cmd.label = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % a selected entity is scanned
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; % the scan is performed for the current entity
            end
            if length(control.directives(d).options) > 2 % an output filename is provided
                cmd.fname = control.directives(d).options{3};
            else
                cmd.fname = 'MMMx_site_scan'; % default file name
            end
            if length(control.directives(d).options) > 3 % residue type selection is made
                cmd.restypes = control.directives(d).options{4};
            else
                cmd.restypes = residue_types; % default residue types
            end
            if length(control.directives(d).options) > 4 % minimum number of rotamers
                cmd.min_rotamers = str2double(control.directives(d).options{5});
            else
                cmd.min_rotamers = min_rotamers; % use default
            end
            if length(control.directives(d).options) > 5 % minimum partition function
                cmd.min_Z = str2double(control.directives(d).options{6});
            else
                cmd.min_Z = min_partition_function; % use default
            end
            if length(control.directives(d).options) > 6 % select only certain chains
                cmd.chains = control.directives(d).options{7};
            else
                cmd.chains = '*'; % all chains
            end
            commands{cmd_poi} = cmd;
        case {'pairlist','pairlist!'}
            cmd_poi = cmd_poi + 1;
            cmd.sitescan = control.directives(d).options{1};
            cmd.entity = control.directives(d).options{2};
            if length(control.directives(d).options) > 2 % output filename is provided
                cmd.fname = control.directives(d).options{3};
            else
                cmd.fname = 'MMMx_pair_list'; % default file name
            end
            if length(control.directives(d).options) > 3 % a minimum mean distance is provided
                cmd.r_min = str2double(control.directives(d).options{4});
            else
                cmd.r_min = min_distance; % use default
            end
            if length(control.directives(d).options) > 4 % a maximum mean distance is provided
                cmd.r_max = str2double(control.directives(d).options{5});
            else
                cmd.r_max = max_distance; % use default
            end
            commands{cmd_poi} = cmd;
        case 'rbreference'
            cmd_poi = cmd_poi + 1;
            cmd.entity = control.directives(d).options{1};
            cmd.r_min = str2double(control.directives(d).options{2});            
            cmd.r_max = str2double(control.directives(d).options{3});            
            cmd.sitescans = cell(1,length(control.directives(d).options)-3);
            for k = 4:length(control.directives(d).options)
                cmd.sitescans{k-3} = control.directives(d).options{k};
            end
            commands{cmd_poi} = cmd;
        case {'hetpairlist','hetpairlist!'}
            cmd_poi = cmd_poi + 1;
            cmd.sitescan1 = control.directives(d).options{1};
            cmd.sitescan2 = control.directives(d).options{2};
            cmd.entity = control.directives(d).options{3};
            if length(control.directives(d).options) > 3 % output filename is provided
                cmd.fname = control.directives(d).options{4};
            else
                cmd.fname = 'MMMx_pair_list'; % default file name
            end
            if length(control.directives(d).options) > 4 % a minimum mean distance is provided
                cmd.r_min = str2double(control.directives(d).options{5});
            else
                cmd.r_min = min_distance; % use default
            end
            if length(control.directives(d).options) > 5 % a maximum mean distance is provided
                cmd.r_max = str2double(control.directives(d).options{6});
            else
                cmd.r_max = max_distance; % use default
            end
            commands{cmd_poi} = cmd;
        case 'distributions'
            cmd_poi = cmd_poi + 1;
            cmd.options.rmin = 10;
            cmd.options.rmax = 150;
            cmd.options.resolution = 0.5;
            cmd.label = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % a selected entity is analyzed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; % distributions are computed for current entity
            end
            if length(control.directives(d).options) > 2 % basis ouput file name is provided
                cmd.fname = control.directives(d).options{3};
            else
                cmd.fname = 'MMMx_distribution'; % default basis output file name
            end
            if length(control.directives(d).options) > 3 % minimum distance is provided
                cmd.options.rmin = str2double(control.directives(d).options{4});
            end
            if length(control.directives(d).options) > 4 % maximum distance is provided
                cmd.options.rmax = str2double(control.directives(d).options{5});
            end
            if length(control.directives(d).options) > 5 % resolution is provided
                cmd.options.resolution = str2double(control.directives(d).options{6});
            end
            [n,narg] = size(control.directives(d).block);
            pair_lists = cell(n,1);
            pl_poi = 0;
            site_pairs = cell(n,2);
            sp_poi = 0;
            % record site lists and site pairs for which distributions
            % should be computed
            for k = 1:n
                if narg > 1 % argument list contains some site pairs
                    if isempty(control.directives(d).block{k,2})
                        pl_poi = pl_poi + 1;
                        pair_lists{pl_poi} = control.directives(d).block{k,1};
                    else
                        sp_poi = sp_poi + 1;
                        site_pairs{sp_poi,1} = control.directives(d).block{k,1};
                        site_pairs{sp_poi,2} = control.directives(d).block{k,2};
                    end
                else % the argument list contains only site lists
                    pl_poi = pl_poi + 1;
                    pair_lists{pl_poi} = control.directives(d).block{k,1};
                end
            end
            cmd.pair_lists = pair_lists(1:pl_poi);
            cmd.site_pairs = site_pairs(1:sp_poi,:);
            commands{cmd_poi} = cmd;            
        case 'trivariate'
            cmd_poi = cmd_poi + 1;
            cmd.label = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % a selected entity is analyzed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; % distributions are computed for current entity
            end
            if length(control.directives(d).options) > 2 % basis ouput file name is provided
                cmd.fname = control.directives(d).options{3};
            else
                cmd.fname = 'MMMx_distribution'; % default basis output file name
            end
            [n,narg] = size(control.directives(d).block);
            site_triples = cell(n,3);
            sp_poi = 0;
            % record site lists and site pairs for which distributions
            % should be computed
            for k = 1:n
                if narg > 1 % argument list contains some site triples
                    sp_poi = sp_poi + 1;
                    site_triples{sp_poi,1} = control.directives(d).block{k,1};
                    site_triples{sp_poi,2} = control.directives(d).block{k,2};
                    site_triples{sp_poi,3} = control.directives(d).block{k,3};
                end
            end
            cmd.site_triples = site_triples(1:sp_poi,:);
            commands{cmd_poi} = cmd;            
        case 'enmpairs'
            cmd_poi = cmd_poi + 1;
            cmd.sitescan = control.directives(d).options{1};
            cmd.entity = control.directives(d).options{2};
            if length(control.directives(d).options) > 2 % output filename is provided
                cmd.fname = control.directives(d).options{3};
            else
                cmd.fname = 'MMMx_enm_pairs'; % default file name
            end
            if length(control.directives(d).options) > 3 % a minimum mean distance is provided
                cmd.r_min = str2double(control.directives(d).options{4});
            else
                cmd.r_min = min_distance; % use default
            end
            if length(control.directives(d).options) > 4 % a maximum mean distance is provided
                cmd.r_max = str2double(control.directives(d).options{5});
            else
                cmd.r_max = max_distance; % use default
            end
            commands{cmd_poi} = cmd;
        case 'rigiflex'
            cmd_poi = cmd_poi + 1;
            cmd.UniProtID = control.directives(d).options{1};
            commands{cmd_poi} = cmd;
        otherwise
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_experiment_design:unknown_directive',...
                'directive %s is unknown',lower(control.directives(d).name));
            record_exception(exceptions{warnings},logfid);
    end
end

ensembles = ensembles(1:ensemble_poi);
commands = commands(1:cmd_poi);

fprintf(logfid,'\n%i commands will be executed on %i structures or ensembles\n\n',cmd_poi,ensemble_poi);

% run the command list

for c = 1:cmd_poi
    cmd = commands{c};
    switch cmd.name
        case 'import'
            ensemble_poi = cmd.ensemble;
            ensemble_descriptor = ensembles{ensemble_poi};
            ensemble_name = ensemble_descriptor.name;
            [entity,exceptions] = get_pdb(cmd.input);
            fprintf(logfid,'\nCurrent ensemble is: %s\n',ensemble_name);
            ensembles = store_ensemble(ensemble_name,entity,ensembles);
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
        case 'getalphafold'
            ensemble_poi = cmd.ensemble;
            ensemble_descriptor = ensembles{ensemble_poi};
            ensemble_name = ensemble_descriptor.name;
            [entity,exceptions] = get_AF(cmd.input);
            if ~isempty(exceptions) && ~isempty(exceptions{1})
                warnings = warnings +1;
                exceptions{warnings} = MException('module_ensembleanalysis:file_does_not_exist',...
                    'AlphaFold prediction for UniProt %s could not be opened',cmd.input);
                record_exception(exceptions{warnings},logfid);
                return
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
        case 'expand'
            ensemble_poi = cmd.ensemble;
            ensemble_descriptor = ensembles{ensemble_poi};
            ensemble_name = ensemble_descriptor.name;
            rba = load(cmd.input);
            entity = get_rba(rba.entity,1); % initialize ensemble entity
            rba_num = length(rba.entity.rba_populations);
            for kr = 2:rba_num
                entity1 = get_rba(rba.entity,kr); % expand next conformer
                fname = sprintf('tmp_MMMx_rba_%i.pdb',kr);
                put_pdb(entity1,fname);
                [entity,exceptions] = get_pdb(fname,[],entity); % add next RBA to ensemble entity
                if ~isempty(exceptions) && ~isempty(exceptions{1})
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_experiment_design:file_does_not_exist',...
                        'PDB file %s could not be opened',fname);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
                delete(fname); % delete temporary file
            end
            ensembles = store_ensemble(ensemble_name,entity,ensembles);
        case 'plot'
            if strcmpi(cmd.extension,'off') || isempty(cmd.fname)
                save_figures = false;
            else
                save_figures = true;
                figure_format = cmd.extension;
                figure_basname = cmd.fname;
            end
        case {'sitescan','sitescan!'}
            c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
            if isempty(c_entity)
                warnings = warnings +1;
                exceptions{warnings} = MException('module_experiment_design:entity_unknown',...
                    'tried a site scan for entity %s, which is unknown',cmd.entity);
                record_exception(exceptions{warnings},logfid);
                return
            end
            clear sitescan_options
            sitescan_options.restypes = cmd.restypes;
            sitescan_options.min_rotamers = cmd.min_rotamers;
            sitescan_options.min_Z = cmd.min_Z;
            sitescan_options.chains = cmd.chains;
            if strcmpi(cmd.name,'sitescan!')
                sitescan_options.ensemble = true;
            else
                sitescan_options.ensemble = false;
            end
            [c_entity,failures] = site_scan_labels(c_entity,cmd.label,cmd.fname,sitescan_options);
            fprintf(logfid,'\nSite scan completed\n');
            fprintf(logfid,'   for %i sites labelling was impossible\n',failures.impossible);
            fprintf(logfid,'   for %i further sites, less than %i rotamer(s) are populated\n',failures.too_few,sitescan_options.min_rotamers);
            fprintf(logfid,'   for %i further sites, partition function was lower than %6.3f\n',failures.too_low_Z,sitescan_options.min_Z);
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            ensembles = store_ensemble(cmd.entity,c_entity,ensembles);
        case 'distributions'
            c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
            if isempty(c_entity)
                warnings = warnings +1;
                exceptions{warnings} = MException('module_experiment_design:entity_unknown',...
                    'tried to compute distance distributions for entity %s, which is unknown',cmd.entity);
                record_exception(exceptions{warnings},logfid);
                return
            end
            if save_figures
                figures = figure_format;
            else
                figures = '';
            end
            labels = split(cmd.label,'|');
            if length(labels) == 1
                labels{2} = labels{1};
            end
            [pairs,~] = size(cmd.site_pairs);
            for p = 1:pairs
               c_entity =  mk_distribution(c_entity,cmd.site_pairs{p,1},labels{1},cmd.site_pairs{p,2},labels{2},...
                   cmd.fname,figures,cmd.fname,cmd.options);
            end
            labels0 = labels;
            for l = 1:length(cmd.pair_lists)
                [pairs,labels] = rd_pair_list(cmd.pair_lists{l});
                if ~strcmp(labels0{1},labels{1}) || ~strcmp(labels0{2},labels{2})
                    fprintf(logfid,'\nIgnoring label specification %s for pair list %s\n',cmd.label,cmd.pair_lists{l});
                    if ~strcmp(labels{1},labels{2})
                        fprintf(logfid,'Labels %s and %s are specified in pair list\n',labels{1},labels{2});
                    else
                        fprintf(logfid,'Label %s is specified in pair list\n',labels{1});
                    end
                end
                for p = 1:length(pairs)
                    c_entity =  mk_distribution(c_entity,pairs(p).address1,labels{1},...
                        pairs(p).address2,labels{2},cmd.fname,figures,cmd.fname,cmd.options);
                end
            end
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            ensembles = store_ensemble(cmd.entity,c_entity,ensembles);
        case 'trivariate'
            c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
            if isempty(c_entity)
                warnings = warnings +1;
                exceptions{warnings} = MException('module_experiment_design:entity_unknown',...
                    'tried to compute trivariate distance distribution for entity %s, which is unknown',cmd.entity);
                record_exception(exceptions{warnings},logfid);
                return
            end
            if save_figures
                figures = figure_format;
            else
                figures = '';
            end
            labels = split(cmd.label,'|');
            if length(labels) == 1
                labels{2} = labels{1};
                labels{3} = labels{1};
            end
            [triples,~] = size(cmd.site_triples);
            for p = 1:triples
               c_entity =  mk_trivariate(c_entity,cmd.site_triples{p,1},labels{1},cmd.site_triples{p,2},labels{2},...
                   cmd.site_triples{p,3},labels{3},cmd.fname,figures,cmd.fname);
            end
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            ensembles = store_ensemble(cmd.entity,c_entity,ensembles);
        case {'pairlist','pairlist!'}
            c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
            if isempty(c_entity)
                warnings = warnings +1;
                exceptions{warnings} = MException('module_experiment_design:entity_unknown',...
                    'tried to generate site pair list for entity %s, which is unknown',cmd.entity);
                record_exception(exceptions{warnings},logfid);
                return
            end
            clear pairlist_options
            pairlist_options.r_min = cmd.r_min;
            pairlist_options.r_max = cmd.r_max;
            if strcmpi(cmd.name,'pairlist!')
                pairlist_options.coupled = false;
            else
                pairlist_options.coupled = true;
            end
            [sites,label] = rd_site_list(cmd.sitescan);
            [c_entity,failures] = pair_list(c_entity,sites,label,cmd.fname,pairlist_options);
            fprintf(logfid,'\nPair list generated\n');
            fprintf(logfid,'   for %i pairs, the distribution is empty in the sensitive range\n',failures.empty);
            fprintf(logfid,'   for %i pairs, mean distance was below %5.1f %s\n',failures.too_short,pairlist_options.r_min,char(197));
            fprintf(logfid,'   for %i pairs, mean distance was above %5.1f %s\n',failures.too_long,pairlist_options.r_max,char(197));
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            ensembles = store_ensemble(cmd.entity,c_entity,ensembles);
        case {'hetpairlist','hetpairlist!'}
            c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
            if isempty(c_entity)
                warnings = warnings +1;
                exceptions{warnings} = MException('module_experiment_design:entity_unknown',...
                    'tried to generate site pair list for entity %s, which is unknown',cmd.entity);
                record_exception(exceptions{warnings},logfid);
                return
            end
            clear pairlist_options
            pairlist_options.r_min = cmd.r_min;
            pairlist_options.r_max = cmd.r_max;
            if strcmpi(cmd.name,'hetpairlist!')
                pairlist_options.coupled = false;
            else
                pairlist_options.coupled = true;
            end
            [sites1,label1] = rd_site_list(cmd.sitescan1);
            [sites2,label2] = rd_site_list(cmd.sitescan2);
            [c_entity,failures] = hetero_pair_list(c_entity,sites1,label1,sites2,label2,cmd.fname,pairlist_options);
            fprintf(logfid,'\nPair list generated\n');
            fprintf(logfid,'   for %i pairs, the distribution is empty in the sensitive range\n',failures.empty);
            fprintf(logfid,'   for %i pairs, mean distance was below %5.1f %s\n',failures.too_short,pairlist_options.r_min,char(197));
            fprintf(logfid,'   for %i pairs, mean distance was above %5.1f %s\n',failures.too_long,pairlist_options.r_max,char(197));
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            ensembles = store_ensemble(cmd.entity,c_entity,ensembles);
        case 'rbreference'
            c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
            if isempty(c_entity)
                warnings = warnings +1;
                exceptions{warnings} = MException('module_experiment_design:entity_unknown',...
                    'tried to generate reference sites for entity %s, which is unknown',cmd.entity);
                record_exception(exceptions{warnings},logfid);
                return
            end
            clear refpoint_options
            refpoint_options.r_min = cmd.r_min;
            refpoint_options.r_max = cmd.r_max;
            all_site_coor_pop(10000).coor = [];
            all_site_coor_pop(10000).pop = [];
            all_site_coor_pop(10000).label = '';
            all_site_coor_pop(10000).address = '';
            poi = 0;
            chains = '';
            for k = 1:length(cmd.sitescans)
                [sites,labels,chains0] = rd_site_list(cmd.sitescans{k});
                if ~strcmp(chains0,'*')
                    chains = strcat(chains,chains0);
                end
                [c_entity,site_coor_pop] = site_geometries(c_entity,sites,labels);
                nsites = length(site_coor_pop);
                all_site_coor_pop(poi+1:poi+nsites) = site_coor_pop;
                poi = poi + nsites;
            end
            all_site_coor_pop = all_site_coor_pop(1:poi);
            [refpoints,failed] = get_refpoints(all_site_coor_pop,refpoint_options);
            fprintf(logfid,'\nReference points in entity %s using site scan lists\n',cmd.entity);
            for k = 1:length(cmd.sitescans)
                fprintf(logfid,'   %s\n',cmd.sitescans{k});
            end
            fprintf(logfid,'\n');
            if failed
                fprintf(logfid,'### No set of three sites found with (rmin,rmax) = (%4.1f,%4.1f) %s ###\n',...
                    refpoint_options.r_min,refpoint_options.r_max,char(197));
            else
                chain_tags = '';
                for k = 1:length(chains)
                    chain_tags = sprintf('%s(%c) ',chain_tags,chains(k));
                end
                fprintf(logfid,'   rigid %s\n',chain_tags);
                for k = 1:3
                    fprintf(logfid,'      %s %s\n',refpoints(k).site,refpoints(k).label);
                end
                fprintf(logfid,'   .rigid\n');
            end
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            ensembles = store_ensemble(cmd.entity,c_entity,ensembles);
        case 'enmpairs'
            c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
            if isempty(c_entity)
                warnings = warnings +1;
                exceptions{warnings} = MException('module_experiment_design:entity_unknown',...
                    'tried to generate ENM pair list for entity %s, which is unknown',cmd.entity);
                record_exception(exceptions{warnings},logfid);
                return
            end
            [sites,label] = rd_site_list(cmd.sitescan);
            labels{1} = label;
            labels{2} = label;
            tic,
            [c_entity,pairs] = score_enm_pairs(c_entity,sites,labels,cmd.fname,cmd.r_min,cmd.r_max);
            toc,
            if ~isfield(pairs,'score')
                fprintf(logfid,'\n### Pair scoring for enm modeling failed ###\n');
            else
                fprintf(logfid,'\n %i site pairs were scored for enm modeling\n',length(pairs));
            end
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            ensembles = store_ensemble(cmd.entity,c_entity,ensembles);
        case 'rigiflex'
            cmd_options.keep_file = true;
            [c_entity,exceptions] = get_AF(cmd.UniProtID,cmd_options);
            if isempty(c_entity)
                warnings = warnings +1;
                exceptions{warnings} = MException('module_experiment_design:no_AlphaFold_prediction',...
                    'AlphaFold prediction does not exist for UniProt ID %s',cmd.UniProtID);
                record_exception(exceptions{warnings},logfid);
                return
            end
            c_entity = domain_partitioning(c_entity);
            ensembles = store_ensemble(cmd.UniProtID,c_entity,ensembles);
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

function entity = mk_distribution(entity,site1,label1,site2,label2,fname,figures,figname,options)

[~,site1] = split_conformer_site(site1);
[~,site2] = split_conformer_site(site2);
tag = sprintf('%s;%s-%s;%s',site1,label1,site2,label2);
filename = sprintf('%s_%s_%s-%s_%s.csv',fname,site1,label1,site2,label2);
figfilename = sprintf('%s_%s_%s-%s_%s.%s',figname,site1,label1,site2,label2,figures);
[r_axis,distribution,entity] = distance_distribution(entity,site1,label1,site2,label2,options);
if ~isempty(r_axis) && ~isempty(distribution) 
    distribution = distribution/sum(distribution);
    distribution = distribution/(r_axis(2)-r_axis(1));
    r_cut = r_axis(distribution > 0.002*max(distribution));
    r_span = max(r_cut) - min(r_cut);
    r_min = min(r_cut) - 0.25*r_span;
    r_max = max(r_cut) + 0.25*r_span;
    h = figure;
    plot(r_axis,distribution,'k','LineWidth',1);
    xlabel(sprintf('distance (%s)',char(197)));
    ylabel(sprintf('probability density (%s^{-1})',char(197)));
    title(tag);
    axis([r_min,r_max,-0.05*max(distribution),1.05*max(distribution)]);
    data = [ r_axis' distribution'];
    writematrix(data,filename);
    if ~isempty(figures)
        saveas(h,figfilename);
    end
end

function entity = mk_trivariate(entity,site1,label1,site2,label2,site3,label3,fname,figures,figname)

[~,site1] = split_conformer_site(site1);
[~,site2] = split_conformer_site(site2);
[~,site3] = split_conformer_site(site3);
filename = sprintf('%s_%s_%s-%s_%s-%s_%s.mat',fname,site1,label1,site2,label2,site3,label3);
[r_axes,trivariate,entity] = trivariate_distribution(entity,site1,label1,site2,label2,site3,label3);
if ~isempty(r_axes) && ~isempty(trivariate) 
    % normalize trivariate distribution
    trivariate = trivariate/sum(sum(sum(trivariate)));
    voxel_volume = (r_axes{1}(2)-r_axes{1}(1))*(r_axes{2}(2)-r_axes{2}(1))*(r_axes{3}(2)-r_axes{3}(1));
    trivariate = trivariate/voxel_volume;
    % save distribution and axes
    r_axis_1 = r_axes{1};
    r_axis_2 = r_axes{2};
    r_axis_3 = r_axes{3};
    save(filename,'r_axis_1','r_axis_2','r_axis_3','trivariate');
    % correlation a,b
    tag = sprintf('%s.%s-%s.%s|%s.%s-%s.%s',site1,label1,site2,label2,site1,label1,site3,label3);
    figfilename = sprintf('%s_%s_%s-%s_%s_corr_%s_%s-%s_%s.%s',figname,site1,label1,site2,label2,site1,label1,site3,label3,figures);
    h = figure;
    bivariate = squeeze(sum(trivariate,3));
    contourf(r_axis_1,r_axis_2,bivariate);
    colormap(hot);
    colormap(flipud(colormap));
    title(tag);
    xlabel(sprintf('distance %s.%s (%s)',site1,site2,char(197)));
    ylabel(sprintf('distance %s.%s (%s)',site1,site3,char(197)));
    if ~isempty(figures)
        saveas(h,figfilename);
    end
    tag = sprintf('%s.%s-%s.%s',site1,label1,site2,label2);
    figfilename = sprintf('%s_%s_%s-%s_%s.%s',figname,site1,label1,site2,label2,figures);
    h = figure;
    univariate = sum(bivariate,2);
    univariate = univariate/sum(univariate);
    univariate = univariate/(r_axis_1(2)-r_axis_1(1));
    plot(r_axis_1,univariate,'k','LineWidth',1);
    title(tag);
    xlabel(sprintf('distance %s.%s (%s)',site1,site2,char(197)));
    ylabel(sprintf('probability density (%s^{-1})',char(197)));
    if ~isempty(figures)
        saveas(h,figfilename);
    end
    figfilename = sprintf('%s_%s_%s-%s_%s.%s',figname,site1,label1,site3,label3,figures);
    h = figure;
    univariate = sum(bivariate,1);
    univariate = univariate/sum(univariate);
    univariate = univariate/(r_axis_2(2)-r_axis_2(1));
    plot(r_axis_2,univariate,'k','LineWidth',1);
    title(tag);
    xlabel(sprintf('distance %s.%s (%s)',site1,site3,char(197)));
    ylabel(sprintf('probability density (%s^{-1})',char(197)));
    if ~isempty(figures)
        saveas(h,figfilename);
    end
    % correlation a,c
    tag = sprintf('%s.%s-%s.%s|%s.%s-%s.%s',site1,label1,site2,label2,site2,label2,site3,label3);
    figfilename = sprintf('%s_%s_%s-%s_%s_corr_%s_%s-%s_%s.%s',figname,site1,label1,site2,label2,site2,label2,site3,label3,figures);
    h = figure;
    bivariate = squeeze(sum(trivariate,2));
    contourf(r_axis_1,r_axis_3,bivariate);
    colormap(hot);
    colormap(flipud(colormap));
    title(tag);
    xlabel(sprintf('distance %s.%s (%s)',site1,site2,char(197)));
    ylabel(sprintf('distance %s.%s (%s)',site2,site3,char(197)));
    if ~isempty(figures)
        saveas(h,figfilename);
    end
    tag = sprintf('%s.%s-%s.%s',site2,label2,site3,label3);
    figfilename = sprintf('%s_%s_%s-%s_%s.%s',figname,site2,label2,site3,label3,figures);
    h = figure;
    univariate = sum(bivariate,1);
    univariate = univariate/sum(univariate);
    univariate = univariate/(r_axis_3(2)-r_axis_3(1));
    plot(r_axis_3,univariate,'k','LineWidth',1);
    title(tag);
    xlabel(sprintf('distance %s.%s (%s)',site2,site3,char(197)));
    ylabel(sprintf('probability density (%s^{-1})',char(197)));
    if ~isempty(figures)
        saveas(h,figfilename);
    end
    % correlation b,c
    tag = sprintf('%s.%s-%s.%s|%s.%s-%s.%s',site1,label1,site3,label3,site2,label2,site3,label3);
    figfilename = sprintf('%s_%s_%s-%s_%s_corr_%s_%s-%s_%s.%s',figname,site1,label1,site3,label3,site2,label2,site3,label3,figures);
    h = figure;
    bivariate = squeeze(sum(trivariate,1));
    contourf(r_axis_2,r_axis_3,bivariate);
    colormap(hot);
    colormap(flipud(colormap));
    title(tag);
    xlabel(sprintf('distance %s.%s (%s)',site1,site3,char(197)));
    ylabel(sprintf('distance %s.%s (%s)',site2,site3,char(197)));
    if ~isempty(figures)
        saveas(h,figfilename);
    end
end

function [c,site] = split_conformer_site(address)
% retrieves conformer number from site address

c = 1;
site = address;
pa = strfind(address,'{');
pe = strfind(address,'}');
if ~isempty(pa) && ~isempty(pe)
    c = str2double(address(pa+1:pe-1));
    if pa > 1
        pre = address(1:pa-1);
    else
        pre = '';
    end
    if pe < length(address)
        past = address(pe+1:end);
    else
        past = '';
    end
    site = [pre past];
end

function [c_entity,site_coor_pop] = site_geometries(c_entity,sites,label)

site_coor_pop(length(sites)).coor = [];
site_coor_pop(length(sites)).pop = [];
site_coor_pop(length(sites)).label = '';
site_coor_pop(length(sites)).address = '';

poi = 0;
for k = 1:length(sites)
    [argsout,c_entity,exceptions] = get_label(c_entity,label,{'positions','populations'},sites(k).address);
    if ~isempty(exceptions{1})
        continue
    end
    poi = poi + 1;
    site_coor_pop(poi).coor = argsout{1}{1};
    site_coor_pop(poi).pop = argsout{2}{1};
    site_coor_pop(poi).label = label;
    site_coor_pop(poi).address = sites(k).address;
end
site_coor_pop = site_coor_pop(1:poi);

function [refpoints,failed] = get_refpoints(all_site_coor_pop,refpoint_options)

S = length(all_site_coor_pop); % total number of sites
max_area = 0;
best_sites = zeros(1,3);
% check all triangles
for s1 = 1:S-2
    coor1 = all_site_coor_pop(s1).coor;
    pop1 = all_site_coor_pop(s1).pop;
    for s2 = s1+1:S-1
        coor2 = all_site_coor_pop(s2).coor;
        pop2 = all_site_coor_pop(s2).pop;
        a = get_r_mean(coor1,pop1,coor2,pop2);
        if a < refpoint_options.r_min || a > refpoint_options.r_max
            continue
        end
        for s3 = s2+1:S
            coor3 = all_site_coor_pop(s3).coor;
            pop3 = all_site_coor_pop(s3).pop;
            b = get_r_mean(coor1,pop1,coor3,pop3);
            if b < refpoint_options.r_min || b > refpoint_options.r_max
                continue
            end
            c = get_r_mean(coor2,pop2,coor3,pop3);
            if c < refpoint_options.r_min || c > refpoint_options.r_max
                continue
            end
            s = (a+b+c)/2; % half the circumference
            area = sqrt(s*(s-a)*(s-b)*(s-c)); % area of the triangle, Heron's theorem
            if area > max_area
                max_area = area;
                best_sites = [s1,s2,s3];
            end
        end
    end
end

refpoints(1).site = '';
refpoints(2).site = '';
refpoints(3).site = '';
refpoints(1).label = '';
refpoints(2).label = '';
refpoints(3).label = '';
if prod(best_sites) == 0
    failed = true;
    return
else
    failed = false;
    for k = 1:3
        refpoints(k).site = all_site_coor_pop(best_sites(k)).address;
        refpoints(k).label = all_site_coor_pop(best_sites(k)).label;
    end
end

function r_mean = get_r_mean(positions1,populations1,positions2,populations2)
% Computes a pair distance distribution as well as the fraction of missing
% population, the distribution is normalized to unity

% normalize population vectors
populations1 = populations1/sum(populations1);
populations2 = populations2/sum(populations2);

[n_rot_1,~] = size(positions1);   % returns number of rotamers at position 1
[n_rot_2,~] = size(positions2);   % -//- at position 2

a2 = repmat(sum(positions1.^2,2),1,n_rot_2);
b2 = repmat(sum(positions2.^2,2),1,n_rot_1).';
pair_dist = sqrt(abs(a2 + b2 - 2*positions1*positions2.'));
weights = populations1*populations2.';
r_mean = sum(sum(pair_dist.*weights));