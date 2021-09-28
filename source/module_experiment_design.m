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
        case {'addpdb','getens','expand','import'}
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
        case 'distributions'
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
            site_lists = cell(n,1);
            sl_poi = 0;
            site_pairs = cell(n,2);
            sp_poi = 0;
            % record site lists and site pairs for which distributions
            % should be computed
            for k = 1:n
                if narg > 1 % argument list contains some site pairs
                    if isempty(control.directives(d).block{k,2})
                        sl_poi = sl_poi + 1;
                        site_lists{sl_poi} = control.directives(d).block{k,1};
                    else
                        sp_poi = sp_poi + 1;
                        site_pairs{sp_poi,1} = control.directives(d).block{k,1};
                        site_pairs{sp_poi,2} = control.directives(d).block{k,2};
                    end
                else % the argument list contains only site lists
                    sl_poi = sl_poi + 1;
                    site_lists{sl_poi} = control.directives(d).block{k,1};
                end
            end
            cmd.site_lists = site_lists(1:sl_poi);
            cmd.site_pairs = site_pairs(1:sp_poi,:);
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
            [c_entity,failures] = site_scan(c_entity,cmd.label,cmd.fname,sitescan_options);
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
    end
end
disp('Aber hallo!');

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


function h = plot_pair_rmsd(pair_rmsd,superimposed)

if ~exist('superimposed','var') || isempty(superimposed)
    superimposed = true;
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
if superimposed
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
legend([h1,h2,h3],'segment length distribution','mean value',sprintf('random coil %5.3f k^{%5.3f}',segments.R0_ee,segments.nu_ee),'Location','southeast');
axis([min(kaxis)-1,max(kaxis)+1,0,1.05*max(segments.max_R2)]);