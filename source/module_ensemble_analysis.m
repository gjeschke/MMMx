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
        case 'acs' % abstract conformer space
            cmd_poi = cmd_poi + 1;
            arg1 = split(control.directives(d).options{1},'.');
            cmd.entities{1} = arg1{1};
            if length(arg1) > 1
                cmd.addresses{1} = arg1{2};
            else
                cmd.addresses{1} = '';
            end
            for k = 2:length(control.directives(d).options)
                argk =split(control.directives(d).options{k},'.');
                cmd.entities{k} = argk{1};
                if length(argk) > 1
                    cmd.addresses{k} = argk{2};
                else
                    cmd.addresses{k} = '';
                end
            end
            commands{cmd_poi} = cmd;
        case {'addpdb','getens','input','get_ped','get_zenodo','get_mmmx'}
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
        case 'zenodo'
            cmd_poi = cmd_poi + 1;
            cmd.input = control.directives(d).options{1};
            commands{cmd_poi} = cmd;
        case 'cluster'
            cmd_poi = cmd_poi + 1;
            cmd.entity = control.directives(d).options{1};
            cmd.new_ensemble = control.directives(d).options{2};
            cmd.size = str2double(control.directives(d).options{3});
            cmd.address = '';
            if length(control.directives(d).options) > 3 % chain and possibly range given
                cmd.address = control.directives(d).options{4};
            end            
            commands{cmd_poi} = cmd;
        case 'transition'
            cmd_poi = cmd_poi + 1;
            cmd.entity1 = control.directives(d).options{1};
            cmd.entity2 = control.directives(d).options{2};
            cmd.superimposed = '';
            cmd.visualize = '';
            if length(control.directives(d).options) >= 3
                cmd.superimposed = control.directives(d).options{3};
            end
            if length(control.directives(d).options) >= 4
                cmd.visualize = control.directives(d).options{4};
            end
            cmd.graphics = control.directives(d).block;
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
        case 'match'
            cmd_poi = cmd_poi + 1;
            cmd.entity1 = control.directives(d).options{1};
            cmd.entity2 = control.directives(d).options{2};
            cmd.address = '';
            cmd.resolved = false;
            if length(control.directives(d).options) > 2 % chain and possibly range given
                cmd.address = control.directives(d).options{3};
            end
            if length(control.directives(d).options) > 3 % separate chain and possibly range given for second ensemble
                cmd.address2 = control.directives(d).options{4};
            else
                cmd.address2 = cmd.address;
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
            cmd.drms = true;
            cmd.maxpop = false;
            cmd.population = false;
            if length(control.directives(d).options) > 1 % a selected entity is analyzed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; 
            end
            if length(control.directives(d).options) > 2 % further directive
                if strcmpi(control.directives(d).options{3},'drms')
                    cmd.drms = true;
                end
                if strcmpi(control.directives(d).options{3},'similarity')
                    cmd.drms = true;
                    cmd.maxpop = true;
                end
                if strcmpi(control.directives(d).options{3},'oriented')
                    cmd.oriented = true;
                    cmd.drms = false;
                end
                if strcmpi(control.directives(d).options{3},'population')
                    cmd.oriented = true;
                    cmd.drms = false;
                    cmd.population = true;
                end
            end
            commands{cmd_poi} = cmd;
        case {'save','archive','put_mmmx'}
            cmd_poi = cmd_poi + 1;
            cmd.outname = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % a selected entity is saved
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
            cmd.options.pair_drms = true;
            cmd.options.superimpose = true;
            cmd.options.sorted = false;
            cmd.options.pair_corr = false;
            cmd.options.compactness = false;
            cmd.options.matlab = false;
            cmd.options.csv = false;
            cmd.options.chain = '';
            cmd.options.range = [];
            cmd.options.compactness_range = [];
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
                        [~,args] = size(control.directives(d).block);
                        if args >= 3
                            range_min = str2double(control.directives(d).block{k,2});
                            range_max = str2double(control.directives(d).block{k,3});
                            if isnan(range_max) && ~isnan(range_min)
                                range_max = range_min;
                                range_min = -range_max;
                            end
                            if ~isnan(range_max) && ~isnan(range_min)
                                cmd.options.compactness_range = [range_min,range_max];
                            end
                        elseif args == 2
                            range_max = str2double(control.directives(d).block{k,3});
                            if ~isnan(range_max)
                                cmd.options.compactness_range = [-range_max,range_max];
                            end
                        end
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
            cmd.options.atoms = 'CA'; % do not use 'all', this creates problems with varying protonatioon of sidegroups
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
                cmd.address = control.directives(d).options{2};
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
        case 'acs'
            entities = cell(1,length(cmd.entities));
            entity_names = '';
            for ent = 1:length(cmd.entities)
                entities{ent} = retrieve_ensemble(cmd.entities{ent},ensembles,logfid);
                entity_names = strcat(entity_names,cmd.entities{ent});
                entity_names = strcat(entity_names,'_');
            end
            entity_names = entity_names(1:end-1);
            [coor,coor_3D,measures] = abstract_conformer_space(entities,cmd.addresses);
            datname = sprintf('acs_%s.csv',entity_names);
            writematrix(coor,datname);
            datname = sprintf('acs_3D_%s.csv',entity_names);
            writematrix(coor_3D,datname);
            fprintf(logfid,'--- Abstract conformer analysis for %s',cmd.entities{1});
            for ent = 2:length(cmd.entities)
                fprintf(logfid,', %s',cmd.entities{ent});
            end
            fprintf(logfid,' ---\n\n');
            fprintf(logfid,'Real space radius of gyration: %5.1f %c\n',measures.Rg,char(197));
            fprintf(logfid,'Abstract conformer space radius of gyration: %5.1f %c\n',measures.Rg_acs,char(197));
            fprintf(logfid,'Ensemble disorder: %6.4f\n',measures.disorder);
            fprintf(logfid,'Ensemble extension: %5.1f %c\n',measures.extension,char(197));
            fprintf(logfid,'ACS dimension is %i for %i conformers\n',measures.dimension,length(measures.populations));
            fprintf(logfid,'ACS embedding error is %6.4f %c\n',measures.error_nD,char(197));
            fprintf(logfid,'ACS 3D embedding error is %5.1f %c\n',measures.error_3D,char(197));
            visualize_acs(coor_3D,measures.errors_3D,measures.assignment,measures.populations,entity_names);
        case 'zenodo'
            args = split(cmd.input,'.');
            fname = args{2};
            k = 2;
            while k < length(args)
                k = k + 1;
                fname = sprintf('%s.%s',fname,args{k});
            end
            query = sprintf('https://zenodo.org/api/records/%s',args{1});
            try
                zenodo_info = webread(query);
            catch
                fprintf(logfid,'ERROR: Access to Zenodo ID %s failed.\n',args{1});
            end            
            query = '';
            k = 0;
            while k < length(zenodo_info.files)
                k = k + 1;
                if strcmpi(zenodo_info.files(k).key,fname)
                    query = zenodo_info.files(k).links.self;
                    break
                end
            end            
            try
                websave(fname,query);
            catch
                fprintf(logfid,'ERROR: Download of file %s from Zenodo ID %s failed.\n',fname,args{1});
            end            
            fprintf(logfid,'\nDownloaded file %s from Zenodo ID %s.\n',fname,args{1});
            % unzip file if required
            [~,~,ext] = fileparts(fname);            
            switch ext
                case '.zip'
                    unzip(fname);
                    fprintf(logfid,'Archive %s was unzipped.\n',fname);
                case '.gz'
                    filenames = gunzip(fname);
                    fprintf(logfid,'Archive %s was unzipped.\n',fname);
                    if length(filenames) == 1
                        [~,~,ext2] = fileparts(filenames{1});
                        if strcmpi(ext2,'.tar')
                            untar(filenames{1});
                            fprintf(logfid,'Archive %s was unzipped.\n',filenames{1});
                        end
                    end
                case '.tar'
                    untar(fname);
                    fprintf(logfid,'Archive %s was unzipped.\n',fname);
            end
        case 'addpdb'
            ensemble_poi = cmd.ensemble;
            ensemble_descriptor = ensembles{ensemble_poi};
            ensemble_name = ensemble_descriptor.name;
            [path,~,ext] = fileparts(cmd.input);
            if isempty(ext)
                cmd.input = strcat(cmd.input,'.pdb');
            end
            added_files = dir(cmd.input); % find all files that match the pattern
            filenames = cell(1,length(added_files));
            for conf = 1:length(added_files)
                filenames{conf} = fullfile(path,added_files(conf).name);
            end
            [entity,exceptions] = entity_from_filelist(filenames);
            if ~isempty(exceptions) && ~isempty(exceptions{1})
                return
            end
            ensemble_descriptor.entity = entity;
            fprintf(logfid,'\nCurrent ensemble is: %s\n',ensemble_name);
            ensembles = store_ensemble(ensemble_name,entity,ensembles);
        case {'getens','input'}
            ensemble_poi = cmd.ensemble;
            ensemble_descriptor = ensembles{ensemble_poi};
            ensemble_name = ensemble_descriptor.name;
            % allow for input of zipped ensembles 
            [~,~,extension] = fileparts(cmd.input);
            if strcmpi(extension,'.zip')
                filenames = unzip(cmd.input);
                for f = 1:length(filenames)
                    [~,~,ext] = fileparts(filenames{f});
                    if strcmpi(ext,'.ens')
                        cmd.input = filenames{f};
                    end
                end
            end
            entity = get_ensemble(cmd.input);
            fprintf(logfid,'\nCurrent ensemble is: %s\n',ensemble_name);
            ensembles = store_ensemble(ensemble_name,entity,ensembles);
        case {'get_ped'}
            ensemble_poi = cmd.ensemble;
            ensemble_descriptor = ensembles{ensemble_poi};
            ensemble_name = ensemble_descriptor.name;
            args = split(cmd.input,'.');
            entity = get_PED(args{1},args{2});
            fprintf(logfid,'\nCurrent ensemble is: %s\n',ensemble_name);
            ensembles = store_ensemble(ensemble_name,entity,ensembles);
        case {'get_zenodo'}
            ensemble_poi = cmd.ensemble;
            ensemble_descriptor = ensembles{ensemble_poi};
            ensemble_name = ensemble_descriptor.name;
            args = split(cmd.input,'.');
            fname = args{2};
            k = 2;
            while k < length(args)
                k = k + 1;
                fname = sprintf('%s.%s',fname,args{k});
            end
            entity = get_zenodo(args{1},fname);
            fprintf(logfid,'\nCurrent ensemble is: %s\n',ensemble_name);
            ensembles = store_ensemble(ensemble_name,entity,ensembles);
        case {'get_mmmx'}
            ensemble_poi = cmd.ensemble;
            ensemble_descriptor = ensembles{ensemble_poi};
            ensemble_name = ensemble_descriptor.name;
            temp_data = load(cmd.input);
            entity = temp_data.entity;
            fprintf(logfid,'\nCurrent ensemble is: %s\n',ensemble_name);
            ensembles = store_ensemble(ensemble_name,entity,ensembles);
            clear temp_data
        case {'cluster'}
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'entity %s cannot be clustered, since entity is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            cluster_options.size = cmd.size;
            if ~isempty(cmd.address)
                [chain,range] = split_chain_range(cmd.address);
                cluster_options.chain = chain;
                cluster_options.range = range;
            end
            [ensemble,info,assignment] = cluster_ensemble(c_entity,cluster_options);
            original = 1:length(assignment);
            [C,~] = size(ensemble);
            fprintf(logfid,'\n--- Ensemble %s clustered to ensemble %s with %i conformers\n\n',cmd.entity,cmd.new_ensemble,C);
            if ~isempty(cmd.address)
                fprintf(logfid,'Clustering with respect to ');
                if ~isempty(chain)
                    fprintf(logfid,'chain %s and ',cluster_options.chain);
                end
                if ~isempty(range)
                    if min(range) > 0
                        fprintf(logfid,'residue range %i-%i',cluster_options.range);
                    else
                        fprintf(logfid,'residues %i,',-cluster_options.range(1));
                        for res = 2:length(cluster_options.range)-1
                            fprintf(logfid,'%i,',-cluster_options.range(res));
                        end
                        fprintf(logfid,'%i',-cluster_options.range(end));
                    end
                end
                fprintf(logfid,'\n');
            end
            fprintf(logfid,'Ensemble entropy reduced from %5.2f to %5.2f\n',info.entropies);
            fprintf(logfid,'Ensemble width changed from %4.1f to %4.1f %c\n',info.widths,char(197));
            fprintf(logfid,'Clustering resolution is %4.1f %c\n',info.resolution,char(197));
            fprintf(logfid,'Similarity of reduced ensemble to original ensemble is %6.4f\n',info.similarity);
            ensemble_name = strcat(cmd.new_ensemble,'.ens');
            ens_fid = fopen(ensemble_name,'wt');
            fprintf(ens_fid,'%% Clustered ensemble %s by MMMx derived from ensemble %s\n\n',cmd.new_ensemble,cmd.entity);
            for clust = 1:C
                fprintf(logfid,'Cluster %i with population %6.4f is represented by conformer %i of original ensemble\n',clust,ensemble(clust,2),ensemble(clust,1));
                members = original(assignment == clust);
                fprintf(logfid,'  included conformers: ');
                for m = 1:length(members)-1
                    fprintf(logfid,'%i,',members(m));
                end
                fprintf(logfid,'%i\n',members(end));
                oname = sprintf('%s_m%i.pdb',cmd.new_ensemble,clust);
                clear save_options
                save_options.order = ensemble(clust,1);
                exceptions = put_pdb(c_entity,oname,save_options);
                fprintf(ens_fid,'%s  %8.6f\n',oname,ensemble(clust,2));
            end
            fclose(ens_fid);
            c_entity = get_ensemble(ensemble_name);
            if strcmp(cmd.entity,'.')
                entity = c_entity;
            end
            ensembles = store_ensemble(cmd.new_ensemble,c_entity,ensembles); 
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
           fprintf(logfid,'\n*** We advise to use the keyword "match" for ensemble comparison ***\n\n');
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
               title(sprintf('\nOverlap between ensembles %s and %s',cmd.entity1,cmd.entity2));
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
               fprintf(logfid,'\nEnsemble density overlap between %s and %s is %6.3f\n',cmd.entity1,cmd.entity2,overlap);
           end
        case 'match'
           entity1 = retrieve_ensemble(cmd.entity1,ensembles,logfid);
           if isempty(entity1)
               warnings = warnings +1;
               exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                   'tried comparison with entity %s, which is unknown',cmd.entity1);
               record_exception(exceptions{warnings},logfid);
               return
           end
           entity2 = retrieve_ensemble(cmd.entity2,ensembles,logfid);
           if isempty(entity2)
               warnings = warnings +1;
               exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                   'tried comparison with entity %s, which is unknown',cmd.entity2);
               record_exception(exceptions{warnings},logfid);
               return
           end
           [chain,range] = split_chain_range(cmd.address);
           [chain2,range2] = split_chain_range(cmd.address2);
           pair_drms = pair_drms_matrix(entity1,chain,range,entity2,chain2,range2);
           mismatch = 0;
           m1 = length(entity1.populations);
           m2 = length(entity2.populations);
           matches = pair_drms(1:m1,m1+1:end);
           selfmatch1 = pair_drms(1:m1,1:m1);
           selfmatch2 = pair_drms(m1+1:end,m1+1:end);
           var11 = kron(entity1.populations,entity1.populations').*selfmatch1.^2;
           var22 = kron(entity2.populations,entity2.populations').*selfmatch2.^2;
           var12 = kron(entity1.populations,entity2.populations').*matches.^2;
           diff = sqrt(sum(sum(var12)) - sum(sum(var11))/2 - sum(sum(var22))/2);
           width1 = sqrt(sum(sum(var11))); 
           width2 = sqrt(sum(sum(var22)));
           similarity = sqrt((width1*width2)/sum(sum(var12)));
           ens1 = cmd.entity1;
           if ~isempty(chain)
               ens1 = sprintf('%s (%s)',ens1,chain);
           end
           if ~isempty(range)
               ens1 = sprintf(ens1,'%s %i-%i',range);
           end
           ens2 = cmd.entity2;
           if ~isempty(chain2)
               ens2 = sprintf('%s (%s)',ens2,chain2);
           end
           if ~isempty(range2)
               ens2 = sprintf(ens2,'%s %i-%i',range2);
           end
           fprintf(logfid,'\nEnsembles %s and %s have a mean variability deviation of %4.1f %c at ensemble widths of %4.1f and %4.1f %c\n',ens1,ens2,diff,char(197),width1,width2,char(197));
           fprintf(logfid,'The similarity measure is %5.3f\n\n',similarity);
           fprintf(logfid,'\n--- Matching conformers in %s by conformers in %s ---\n\n',cmd.entity1,cmd.entity2);
           if m1 <= m2
               mismatches = zeros(m1,1);
               best = zeros(m1,1);
               for k = 1:m1
                   [colmin,lines] = min(matches);
                   [~,col] = min(colmin);
                   best(lines(col)) = col;
                   mismatches(lines(col)) = matches(lines(col),col);
                   mismatch = mismatch + entity1.populations(lines(col))*mismatches(lines(col))^2;
                   matches(lines(col),:) = 1e12;
                   matches(:,col) = 1e12;
               end
               fprintf(logfid,'Mean mismatch is %4.1f %c\n',sqrt(mismatch/sum(entity1.populations)),char(197));
               for k = 1:m1
                   fprintf(logfid,'%s.%i is matched by %s.%i with DRMSD of %4.1f %c\n',cmd.entity1,k,cmd.entity2,best(k),mismatches(k),char(197));
               end
           end
           if m1 > m2
               fprintf(logfid,'\nSmaller ensemble is %s, reverse matching\n\n',cmd.entity2);
               mismatches = zeros(m2,1);
               best = zeros(m2,1);
               for k = 1:m2
                   [colmin,lines] = min(matches);
                   [~,col] = min(colmin);
                   best(col) = lines(col);
                   mismatches(col) = matches(lines(col),col);
                   mismatch = mismatch + entity2.populations(col)*mismatches(col)^2;
                   matches(lines(col),:) = 1e12;
                   matches(:,col) = 1e12;
               end
               fprintf(logfid,'Mean mismatch is %4.1f %c\n',sqrt(mismatch/sum(entity1.populations)),char(197));
               for k = 1:m2
                   fprintf(logfid,'%s.%i is matched by %s.%i with DRMSD of %4.1f %c\n',cmd.entity1,best(k),cmd.entity2,k,mismatches(k),char(197));
               end
           end
        case 'transition'
           clear options
           args = split(cmd.entity1,'.');
           cmd.entity1 = args{1};
           if length(args)>1
               options.chain1 = args{2};
           else
               options.chain1 = '';
           end
           args = split(cmd.entity2,'.');
           cmd.entity2 = args{1};
           if length(args)>1
               options.chain2 = args{2};
           else
               options.chain1 = '';
           end
           entity1 = retrieve_ensemble(cmd.entity1,ensembles,logfid);
           if isempty(entity1)
               warnings = warnings +1;
               exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                   'tried comparison with entity %s, which is unknown',cmd.entity1);
               record_exception(exceptions{warnings},logfid);
               return
           end
           entity2 = retrieve_ensemble(cmd.entity2,ensembles,logfid);
           if isempty(entity2)
               warnings = warnings +1;
               exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                   'tried to comparison with entity %s, which is unknown',cmd.entity2);
               record_exception(exceptions{warnings},logfid);
               return
           end
           options.visualize = cmd.visualize;
           options.fname1 = cmd.entity1;
           options.fname2 = cmd.entity2;
           options.superimposed = cmd.superimposed;
           options.graphics = cmd.graphics;
           if save_figures
               options.figname = sprintf('transition_%s_%s.%s',cmd.entity1,cmd.entity2,figure_format);
           end
           assignment_name = sprintf('transition_%s_%s_assignment.mat',cmd.entity1,cmd.entity2);
           clusters = cluster_transition(entity1,entity2,options);
           save(assignment_name,'clusters');
           fprintf(logfid,'\n--- Transition analysis ---\n\n');
           fprintf(logfid,'Ensemble 1 (%s) has %i conformers.\n',cmd.entity1,clusters.C1);
           fprintf(logfid,'Ensemble 2 (%s) has %i conformers.\n',cmd.entity2,clusters.C2);
           fprintf(logfid,'Similarity of the two ensembles is %5.3f.\n',clusters.similarity);
           fprintf(logfid,'%i clusters were generated\n',clusters.nc);
           fprintf(logfid,'%i pure ensemble 1 clusters\n',sum(clusters.type == 1));
           fprintf(logfid,'%i pure ensemble 2 clusters\n',sum(clusters.type == 2));
           fprintf(logfid,'%i mixed clusters\n',sum(clusters.type == 0));
           [sorted,sorting] = sort(clusters.type);
           pop_dp = 0;
           pop_cs_initial = 0;
           pop_cs_final = 0;
           pop_if = 0;
           pop1 = entity1.populations;
           pop2 = entity2.populations;
           for clust = 1:clusters.nc
               index = sorting(clust);
               mem1 = clusters.members{index,1};
               c1 = length(mem1);
               cpop1 = 0;
               for kc1 = 1:c1
                   cpop1 = cpop1 + pop1(mem1(kc1));
               end
               mem2 = clusters.members{index,2};
               c2 = length(mem2);
               cpop2 = 0;
               for kc2 = 1:c2
                   cpop2 = cpop2 + pop2(mem2(kc2));
               end
               switch sorted(clust)
                   case 0
                       ctype = 'mixed';
                       pop_cs_initial = pop_cs_initial + cpop1;
                       pop_cs_final = pop_cs_final + cpop2;
                   case 1
                       ctype = 'pure ensemble 1';
                       pop_dp = pop_dp + cpop1;
                   case 2
                       ctype = 'pure ensemble 2';
                       pop_if = pop_if + cpop2;
               end
               fprintf(logfid,'\nCluster %i is of %s type with Rg = %4.1f %c and on similarity scale %5.3f\n',clust,ctype,clusters.Rg(index),char(197),clusters.scale(index));
               switch sorted(clust)
                   case 0
                       fprintf(logfid,'%i members, %i from ensemble 1, %i from ensemble 2\n',c1+c2,c1,c2);
                   case 1
                       fprintf(logfid,'%i members from ensemble 1\n',c1);
                   case 2
                       fprintf(logfid,'%i members from ensemble 2\n',21);
               end
           end
           fprintf(logfid,'\nPopulation of initial state that is deselected is %5.1f\n',100*pop_dp);
           fprintf(logfid,'Population of initial state that is retained is %5.1f\n',100*pop_cs_initial);
           fprintf(logfid,'\nPopulation of final state from conformational selection is %5.1f\n',100*pop_cs_final);
           fprintf(logfid,'Population of final state from induced fit is %5.1f\n',100*pop_if);
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
            if cmd.population
                pop = c_entity.populations;
                [cluster_pop,ordering] = sort(pop,'descend');
                cluster_sizes = ones(1,length(ordering));
            else
                if cmd.drms
                    [pair_rmsd,pop,Rg,exceptions0] = pair_drms_matrix(c_entity);
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
                if cmd.maxpop
                    [pair_rmsd,Rg,ordering] = similarity_sorting(pair_rmsd,Rg);
                    cluster_sizes = ones(1,length(ordering));
                    cluster_pop = pop(ordering);
                    figure(777);
                    plot(1:length(Rg),Rg,'.','MarkerSize',10);
                    xlabel('Conformer number');
                    ylabel(sprintf('Rg (%c)',char(197)));
                else
                    [pair_rmsd,ordering,cluster_assignment,cluster_sizes,cluster_pop] = cluster_sorting(pair_rmsd,pop);
                    D = dunn_index(pair_rmsd,cluster_assignment);
                    fprintf(logfid,'\nCluster assignment has a Dunn index of %5.3f\n',D);
                end
            end
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
            elseif ~cmd.population
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
                        'tried to save entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            save_options.pop = true;
            save_options.order = 1:length(c_entity.populations);
            exceptions0 = put_pdb(c_entity,cmd.outname,save_options);
            for exci = 1:length(exceptions0)
                if ~isempty(exceptions0{exci})
                    warnings = warnings +1;
                    exceptions{warnings} = exceptions0{exci};
                    record_exception(exceptions{warnings},logfid);
                end
            end
            fprintf(logfid,'\nEnsemble %s saved to single PDB file %s\n',cmd.entity,cmd.outname);
            [pname,basname] = fileparts(cmd.outname);
            wname = fullfile(pname, sprintf('%s_weights.tsv',basname));
            fid = fopen(wname,'wt');
            for conf = 1:length(c_entity.populations)
                fprintf(fid,'%i\t%8.6f\n',conf,c_entity.populations(conf));
            end
            fclose(fid);
        case 'archive'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'tried to archive entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            [~,basname,~] = fileparts(cmd.outname);
            filenames = cell(1,1+length(c_entity.populations));
            filenames{1} = strcat(basname,'.ens');
            ens_fid = fopen(filenames{1},'wt');
            fprintf(ens_fid,'%% Ensemble %s (%s)\n',cmd.entity,c_entity.name);
            for conf = 1:length(c_entity.populations)
                filenames{1+conf} = sprintf('%s_m%i.pdb',basname,conf);
                clear save_options
                save_options.order = conf;
                exceptions = put_pdb(c_entity,filenames{1+conf},save_options);
                fprintf(ens_fid,'%s  %8.6f\n',filenames{1+conf},c_entity.populations(conf));
            end
            fclose(ens_fid);
            zip(strcat(basname,'.zip'),filenames);
            fprintf(logfid,'\nEnsemble %s was archived as file %s.zip\n',cmd.entity,basname);
        case 'put_mmmx'
            if strcmp(cmd.entity,'.')
                c_entity = entity;
            else
                c_entity = retrieve_ensemble(cmd.entity,ensembles,logfid);
                if isempty(c_entity)
                    warnings = warnings +1;
                    exceptions{warnings} = MException('module_ensembleanalysis:entity_unknown',...
                        'tried to save entity %s, which is unknown',cmd.entity);
                    record_exception(exceptions{warnings},logfid);
                    return
                end
            end
            entity0 = entity;
            entity = c_entity;
            [pathname,basname,~] = fileparts(cmd.outname);
            outname = fullfile(pathname,strcat(basname,'.mat'));
            save(outname,'entity');
            entity = entity0;
            fprintf(logfid,'\nEnsemble %s saved in MMMx internal format as %s\n',cmd.entity,cmd.outname);
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
            
            
            fprintf(logfid,'\n--- Analysis of ensemble %s ---\n',cmd.entity);
            if cmd.options.chain_mode
                fprintf(logfid,'   Restricted to chain %s and residue range %i-%i\n',cmd.options.chain,cmd.options.range);
            end
            fprintf(logfid,'\n');
            [measures,correlations] = analyze_ensemble(backbones,pop,cmd.options);
            
            fprintf(logfid,'Ensemble Shannon entropy (decadic logarithm): %4.2f\n',measures.entropy);
            parts = fieldnames(measures);
            for p = 1:length(parts)
                part = parts{p};
                if isstrprop(part(1),'upper') || strcmp(part,'all')
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
                        fprintf(logfid,'\nEnsemble width and density for %s (%s)\n',name,part);
                        fprintf(logfid,'   width = %4.1f %s; density %4.1f %s\n',...
                            measures.(part).width,char(197),measures.(part).density,char(197));
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
                            figname = sprintf('residue_pair_CI95_%s_%s.%s',basname,part,figure_format);
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
                        if ~isempty(cmd.options.compactness_range)
                            compactness_range = cmd.options.compactness_range;
                        else
                            compactness_range = 1.05*max([abs(min(min(Rdev))),max(max(Rdev))])*[-1,1];
                        end
                        h = plot_compactness(Rdev,'Segment length deviation',compactness_range,offset);
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
            property_cube(c_entity,cmd.outname,cmd.property,cmd.address,cmd.resolution,cmd.pH,cmd.I);
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
            [c_entity,~,exceptions0] = superimpose_ensemble(c_entity,selected,t_entity,cmd.central,cmd.options);
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
            fprintf(logfid,'\nMean asphericity for ensemble %s: %5.3f\n',cmd.entity,c_entity.asphericity);
            pop = c_entity.populations;
            h = plot_asphericity(c_entity.rg_c,c_entity.asphericity_c,pop);
            ha = h.CurrentAxes;
            ha.Title.String = sprintf('Mean asphericity %5.3f',c_entity.asphericity);
            ha.XLabel.String = sprintf('R_g (%c)',char(197));
            ha.YLabel.String = 'Asphericity';
            if save_figures
                figname = sprintf('Asphericity_%s.%s',cmd.entity,figure_format);
                saveas(h,figname);
            end
            asphericity_matrix = [c_entity.populations c_entity.rg_c c_entity.asphericity_c];
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
if ~found
    ensemble_descriptor.name = name;
    ensemble_descriptor.entity = entity;
    ensembles{length(ensembles)+1} = ensemble_descriptor;
end

function record_exception(exception,logfid)

fprintf(logfid,'### ensembleanalysis exception: %s ###\n',exception.message);

function [chain,range] = split_chain_range(address)

chain = '';
range = [];
if isempty(address)
    return
end

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
if contains(address,',')
    residues = split(address,',');
    range = zeros(1,length(residues));
    for res = 1:length(residues)
        range(res) = -str2double(residues{res});
    end
else
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
    image(2*pair_corr,'CDataMapping','scaled');
    title('C\alpha-C\alpha distance 95% confidence interval');
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

function h = plot_asphericity(Rg,asph,pop)

h = figure;
min_Rg = floor(min(Rg));
max_Rg = ceil(max(Rg));
Rg_ax = min_Rg:max_Rg; % 1 Angstroem resolution
asph_ax = 0:0.01:ceil(100*max(asph))/100; % 0.01 resolution for asphericity
pop = pop/max(pop);

pd = zeros(length(asph_ax),length(Rg_ax));
dec_Rg = 1/(sqrt(2)*(Rg_ax(2)-Rg_ax(1)));
dec_asph = 1/(sqrt(2)*(asph_ax(2)-asph_ax(1)));
for k = 1:length(pop)
    xpsf = exp(-(dec_Rg*(Rg_ax-Rg(k))).^2); % point-spread function along x
    ypsf = exp(-(dec_asph*(asph_ax-asph(k))).^2); % point-spread function along x
    pd = pd + pop(k)*ypsf.'*xpsf;
end
pd = pd/max(max(pd));
levels = 0:0.025:1;
contourf(Rg_ax,asph_ax,pd,levels,'edgecolor','none');
c = colorbar;
c.Label.String = 'Relative probability density';
c.FontSize = 12;
set(gca,'FontSize',12);
xlabel(sprintf('R_g (%c)',char(197)));
ylabel('asphericity');


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