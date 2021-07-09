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
%               .name           'yasara', unused
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
% input         file name of input PDB file
% save [fn]     specify file name fn for output, defaults to mmmx_yasara
% 

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% initialize empty output
exceptions = {[]};
failed = false;
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
output_name = false;

% read restraints
for d = 1:length(control.directives)
    switch lower(control.directives(d).name)
        case 'input'
            inname = control.directives(d).options{1};
            fprintf(logfid,'Structure from PDB file(s) %s will be optimized instead of input entity\n',inname);
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

all_files = dir(inname); % find all files that match the pattern
if isempty(all_files)
    file_list{1} = inname;
else
    file_list = cell(1,length(all_files));
    for c = 1:length(all_files)
        file_list{c} = all_files(c).name;
    end
end

for c = 1:length(file_list)
    inname = file_list{c};
    opt.fname = sprintf('%s_m%i',fname,c);
    temporary_input = false;
    if ~isempty(inname)
        poi = strfind(inname,'.pdb');
        if ~isempty(poi)
            inname = inname(1:poi-1);
        end
        if ~output_name
            opt.fname = strcat(inname,'_yasara');
        end
    else
        if ~exist('entity','var') || isempty(entity)
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_yasara:no_input',...
                    'Neither input entity nor input file are defined.');
            fprintf(logfid,'Errror in Yasara optimizer: Neither input entity nor input file are defined\n');
            failed = true;
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

    % read Yasara result into current entity, if output is requested
    if nargout > 0 && isempty(exceptions{1})
        if ~contains(fname,'.pdb')
            fname = strcat(fname,'.pdb');
        end
        if ~exist(fname,'file')
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_yasara:no_result',...
                'ERROR: Yasara computation failed.');
            entity = [];
        else
            entity = get_pdb(fname);
        end
    end
end