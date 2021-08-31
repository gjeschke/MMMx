function [entity,exceptions,failed] = module_ensemble_analysis(control,logfid)
%
% MODULE_ENSEMBLE_ANALYSIS    Analyses an ensemble
%
%   [entity,exceptions] = MODULE_ENSEMBLE_ANALYSIS(control,logfid)
%   Given a conformer (entity), coordinate transformations or sidechain
%   modifications are applied
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
% figures       format for saving figures, 'off' switches figure saving
%               off, this is the default
% flexibility   residue-specific flexibility by Ramachandran angle
%               distribution
% order         local order of residues with respect to the whole structure 
% getens        get ensemble file
% gettraj       get ensemble from trajectory, requires mdtraj (actually
%               mdconvert.exe on the Matlab path)
% order         order ensemble
% paircorr      pair correlation matrix
% segments      distribution of segment-wise root-mean square end-to-end
%               distances and deviation from random-coil fit
% subsample     reduce ensemble (from a trajectory) by subsampling
% superimpose   superimpose conformers
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
        case {'addpdb','getens'}
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
        case 'symmetry'
            cmd_poi = cmd_poi + 1;
            if ~isempty(control.directives(d).options) && ~isempty(control.directives(d).options{1})
                cmd.mode = control.directives(d).options{1}; % superposition mode [backbone|CA|C4'|all]
            else
                cmd.mode = 'backbone'; % by default, the backbones are superimposed 
            end
            if length(control.directives(d).options) > 1 % a selected entity transformed
                cmd.entity = control.directives(d).options{2};
            else
                cmd.entity = '.'; % symmetry transform is performed for current entity
            end
            [n,~] = size(control.directives(d).block);
            cmd.Cn = n; % order of symmetry axis
            cmd.parts = cell(1,n);
            % store selections in the n symmetry-equivalent parts
            for k = 1:n
                cmd.parts{k} = control.directives(d).block{k,1};
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
                entity.populations = ones(1,length(added_files))/length(added_files);
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
            end
            ensemble_descriptor.entity = entity;
            fprintf(logfid,'\nCurrent ensemble is: %s\n',ensemble_name);
            ensembles = store_ensemble(ensemble_name,entity,ensembles);
        case 'getens'
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
    end
end

function entity = retrieve_ensemble(name,ensembles,logfid)

ensemble = [];
e = 0;
while isempty(ensemble) && e < length(ensembles)
    e = e + 1;
    ensemble_descriptor = ensembles{e};
    if strcmp(name,ensemble_descriptor.name)
        entity = ensemble_descriptor.entity;
    end
end
if isempty(ensemble)
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

