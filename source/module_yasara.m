function [entity,exceptions,failed] = module_yasara(control,logfid,entity)
%
% MODULE_YASARA Performs Yasara computations on an entity
%
%   [exceptions,entity,failed] = MODULE_YASARA(control_file)
%   processes Flex controls and restraints and runs the appropriate Yasara
%   function
%
% INPUT
% control       control structure with fields
%               .name           'Yasara_Refine', unused
%               .options        options struct
%               .directives     array of struct with directives
%               .entity_output  number of the entity to be output
% logfid        file handle for log file
% entity        input entity, can be empty if an input file is specified
%
% OUTPUT
% entity        output entity, empty in NOFIT mode
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
% failed        flag indicating whether the module failed
%
% SUPPORTED DIRECTIVES
%
% input         file name of input PDB file(s)
% addpdb        file name of input PDB file(s), synonym to input
% ensemble      file name of an MMMx ensemble list, overrides input/addpdb
% repack        repack sidegroups by SCWRL4 before running Yasara, may help
%               to prevent Yasara crashes
% console       run Yasara in console mode
% save [fn]     specify file name fn for output, defaults to mmmx_yasara
% 

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% initialize empty output
exceptions = {[]};
failed = false;
repack = false;
warnings = 0;

if ~exist('entity','var')
    entity = [];
end
% set defaults

opt.console = false; % Yasara console mode
if length(control.options) >= 1
    opt.max_time = str2double(control.options{1});
else
    opt.max_time = 1;
end

fname = 'MMMx_yasara';
inname = '';
initial_ensemble = '';
output_name = false;

% read restraints
for d = 1:length(control.directives)
    switch lower(control.directives(d).name)
        case {'input','getpdb'}
            inname = control.directives(d).options{1};
            fprintf(logfid,'Structure from PDB file(s) %s will be optimized\n',inname);
        case 'addpdb'
            inname = control.directives(d).options{1};
            fprintf(logfid,'Structure from PDB file(s) %s will be optimized\n',inname);
        case {'getens','ensemble'}
            initial_ensemble = control.directives(d).options{1};
            fprintf(logfid,'Ensemble %s will be optimized\n',initial_ensemble);
        case 'repack'
            repack = true;
        case 'save'
            fname = control.directives(d).options{1};
            % remove extension if any
            [pathname,filename,~] = fileparts(fname);
            fname = fullfile(pathname,filename);
            output_name = true;
        case 'console'
            opt.console = true;
        otherwise
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_yasara:unknown_directive',...
                'directive %s is unknown',lower(control.directives(d).name));
            record_exception(exceptions{warnings},logfid);

    end
end

pop = [];
if ~isempty(initial_ensemble)
    % allow for input of zipped ensembles
    [~,~,extension] = fileparts(initial_ensemble);
    if strcmpi(extension,'.zip')
        filenames = unzip(initial_ensemble);
        for f = 1:length(filenames)
            [~,~,ext] = fileparts(filenames{f});
            if strcmpi(ext,'.ens')
                initial_ensemble = filenames{f};
            end
        end
    end
    [all_files,pop,exceptions0] = rd_ensemble_definition(initial_ensemble);
    if ~isempty(exceptions0) && ~isempty(exceptions{1})
        warnings = warnings + 1;
        exceptions{warnings} = exceptions0{1};
        record_exception(exceptions{warnings},logfid);
        return
    end
end


if isempty(pop)
    all_files = dir(inname); % find all files that match the pattern
end
if isempty(all_files)
    file_list{1} = inname;
else
    file_list = cell(1,length(all_files));
    for c = 1:length(all_files)
        file_list{c} = all_files(c).name;
    end
end

if ~isempty(pop) % we are refining an ensemble
    efid = fopen(sprintf('%s.ens',fname),'wt');
    fprintf(efid,'%% MMMx ensemble %s refined by Yasara\n',fname);
end
one_failed = false;
for c = 1:length(file_list)
    inname = file_list{c};
    fprintf(logfid,'\nOptimizing %s.\n',inname);
    opt.fname = sprintf('%s_m%i',fname,c);
    temporary_input = false;
    if repack
        entity = get_pdb(inname);
        entity = modify_sidechains(entity,true);
        now_string = datestr(now,'HH-MM-SS');
        inname = sprintf('to_be_optimized_%s_%i.pdb',now_string,round(10000*rand));
        put_pdb(entity,inname);
        temporary_input = true;
    end
    if ~isempty(inname)
        poi = strfind(inname,'.pdb');
        if ~isempty(poi)
            inname1 = inname(1:poi-1);
        end
        if ~output_name
            opt.fname = strcat(inname1,'_yasara');
        end
    else
        if ~exist('entity','var') || isempty(entity)
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_yasara:no_input',...
                    'Neither input entity nor input file are defined.');
            fprintf(logfid,'Errror in Yasara optimizer: Neither input entity nor input file are defined\n');
            failed = true;
            if ~isempty(pop)
                fclose(efid);
            end
            return
        else
            % make a file name that is hopefully unique
            now_string = datestr(now,'HH-MM-SS');
            inname = sprintf('to_be_optimized_%s_%i.pdb',now_string,round(10000*rand));
            put_pdb(entity,inname);
            temporary_input = true;
        end
    end

    exceptions = optimize_by_yasara(logfid,inname,opt);

    if temporary_input
        if exist(inname,'file')
            delete(inname);
        end
    end

    fprintf(logfid,'Saving as %s.pdb\n',opt.fname);

    % read Yasara result into current entity, if output is requested
    if nargout > 0 && isempty(exceptions{1})
        if ~contains(opt.fname,'.pdb')
            fname1 = strcat(opt.fname,'.pdb');
        else
            fname1 = opt.fname;
        end
        if ~exist(fname1,'file')
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_yasara:no_result',...
                'ERROR: Yasara computation failed.');
            one_failed = true;
            entity = [];
        else
            entity = get_pdb(fname1);
            if ~isempty(pop)
                fprintf(efid,'%s%10.6f\n',fname1,pop(c)); 
            end
        end
    end
end
if ~isempty(pop)
    fclose(efid);
end
if one_failed
    fprintf(logfid,'### Output ensemble invalid, because at least one conformer could not be optimized ###\n');
end