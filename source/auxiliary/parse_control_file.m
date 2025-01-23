function [controller,exceptions] = parse_control_file(fname)
%
% PARSE_CONTROL_FILE Parses an MMMx control file for MMModel
%
%   [controller,exceptions] = PARSE_CONTROL_FILE(fname)
%   Returns a controller structure for MMModel, containing tasks with
%   associated restraints and run options
%
% INPUT
% fname     name of the control file, '.mxc' is appended if there is no
%           extension and the file without extension does not exist
%
% OUTPUT
% controller   array (1,number_of_tasks) of control structures
%              each task has two standard subfields:
%              .name      task name
%              .options   run options (can be empty), cell array of strings
%
%              further fields correspond to restraints and directives in a
%              task block
%              the following restraint types are preprocessed:
%              ddr        distance distribution restraints
%              ddr_sym    symmetry-related distance distribution restraints
%              pre        paramagnetic relaxation enhancement restraints
%              saxs       small-angle x-ray scattering restraint
%              sans       small-amgle neutron scattering restraint
%              croslink   crosslink restraints
%              depth      immersion depth restraints
%              alpha      alpha-helix propensity
%              beta       beta-strand propensity
%              polypro    proline-II helix propensity
%              cis        cis-peptide propensity
%              rigid      rigid-body definition
%              
%              restraints that are not listed here are treated as
%              directives (stored unprocessed) and must be processed by
%              modules
% exceptions   error messages and warnings, cell array of MException
%              objects cell containing an empty array, if no exception
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output

controller = [];
exceptions = cell(1,10000);
warnings = 0; % counter for warnings

fid = fopen(fname,'rt');

% treat case when file cannot be opened
if fid == -1
    if ~contains(fname,'.') % there was no extension
        fname = strcat(fname,'.mcx');
        fid = fopen(fname,'rt');
    end
end
if fid == -1
    exceptions = {MException('parse_control_file:cannot_be_opened', 'Control file %s cannot be opened', fname)};
    return
end

% initialize level hierarchy and task sequence
level = 0;
indent = 0;
task = 0;

% read first line
tline = fgetl(fid);
while ischar(tline)
    % remove comments and skip over pure comment line
    comment = strfind(tline,'%');
    if ~isempty(comment)
        tline = tline(1:comment-1);
    end
    if isempty(tline)
        tline = fgetl(fid);
        continue
    end
%     fprintf(1,'%s\n',tline);
    % remove trailing blanks
    tline = deblank(tline);
    % convert 
    % determine new indent and current hierarchy level
    line_length = length(tline);
    tline = strtrim(tline);
    new_indent = line_length - length(tline);
    % skip over empty line
    if isempty(tline)
        tline = fgetl(fid);
        continue
    end
    % retrieve arguments
    arguments = textscan(tline,'%s');
    arguments = arguments{1};
    arg1 = strtrim(char(arguments{1}));
    argument_pointer = 2;
    % allow for block key syntax with colon
    if arg1(1) == ':'
        explicit_block_key = true;
        arg1 = arg1(2:end);
    else
        explicit_block_key = false;
    end
    % determine line type and check level
    switch arg1(1)
        case {'#','!'} % keyword or module call on top level
            if level ~= 0 % we are not on top level, return with error
                exceptions = {MException('parse_control_file:wrong_level_keyword',...
                    'Keyword line appears inside module block: %s', tline)};
                fclose(fid);
                return
            end
            % remove space between #/! and keyword/module name if required
            if length(arg1) == 1
                arg1 = strcat(arg1,strtrim(char(arguments{2})));
                argument_pointer = argument_pointer + 1;
            end
            task = task + 1;
            directive = 0;
            taskname = arg1(2:end);
            controller(task).name = taskname; %#ok<AGROW>
            % store options
            controller(task).options = ...
                cell(1,length(arguments)- argument_pointer +1); %#ok<AGROW>
            for karg =  argument_pointer:length(arguments)
                controller(task).options{karg-argument_pointer+1}...
                    = strtrim(char(arguments{karg}));
            end
        case '.' % block closes
            if level == 0 % we are already on top level, return with error
                exceptions = {MException('parse_control_file:wrong_level_close',...
                    'Attempt to close block on top level: %s', tline)};
                fclose(fid);
                return
            end
            % go up one level
            level = level - 1;
%             fprintf(1,'Gone down one level to %i\n',level)
            % remove module name if required
            if length(arg1) == 1
                argument_pointer = argument_pointer + 1;
            end
            % if a module block closed, store entity numbers for returned
            % entities
            if level == 0
                controller(task).entity_output = ...
                    zeros(1,length(arguments)- argument_pointer +1); %#ok<AGROW>
                for karg =  argument_pointer:length(arguments)
                    entity_number = str2double(char(arguments{karg}));
                    if isnan(entity_number) || abs(entity_number - round(entity_number)) > 0
                        warnings = warnings + 1;
                        exceptions{warnings} = {MException('parse_control_file:non_integer_entity_number',...
                    'Noninteger entity number encountered: %s', char(arguments{karg}))};
                    end
                    controller(task).entity_output(karg-argument_pointer+1)...
                        = entity_number;
                end
            end
        otherwise % this must be level > 0
            % the following led to problems with editors using tabs
            % block keys must be closed with .statements now
            %             if new_indent < indent 
            %                 level = level - 1;
            %             end
            if new_indent > indent && level < 2
                level = level + 1;
%                 fprintf(1,'Gone up one level at indent %i to %i\n',new_indent,level);
            end
            switch level
                case 1
                    directive = directive + 1;
                    argline = 0;
                    controller(task).directives(directive).name = arg1; %#ok<AGROW>
                    % store options
                    controller(task).directives(directive).options = ...
                        cell(1,length(arguments)- argument_pointer +1); %#ok<AGROW>
                    controller(task).directives(directive).modifiers = ...
                        cell(1,length(arguments)- argument_pointer +1); %#ok<AGROW>
                    kopt = 0;
                    kmod = 0;
                    for karg =  argument_pointer:length(arguments)
                        curr_argument = strtrim(char(arguments{karg}));
                        if curr_argument(1) ~= ':'
                            kopt = kopt + 1;
                            controller(task).directives(directive).options{kopt}...
                                = curr_argument;
                        else
                            kmod = kmod + 1;
                            controller(task).directives(directive).modifiers{kmod}...
                                = curr_argument(2:end);
                        end
                    end
                    controller(task).directives(directive).options = ...
                        controller(task).directives(directive).options(1:kopt); %#ok<AGROW>
                    controller(task).directives(directive).modifiers = ...
                        controller(task).directives(directive).modifiers(1:kmod); %#ok<AGROW>
                    controller(task).directives(directive).block = cell(1,1);  %#ok<AGROW>
                case 2
                    argline = argline+1;
                    for karg =  1:length(arguments)
                        controller(task).directives(directive).block{argline,karg}...
                            = strtrim(char(arguments{karg})); %#ok<AGROW>
                    end
            end
    end
    
    indent = new_indent;
    if explicit_block_key
        level = level + 1;
    end
    tline = fgetl(fid);
end


% close control file
fclose(fid);

exceptions = exceptions(1:warnings);
