function [exceptions,entity] = MMModel(control_file)
%
% MMModel Run modelling pipeline
%
%   [exceptions,entity] = MMModel(control_file)
%   Parses control file for modelling directives and arguments and runs the
%   specified modelling pipeline
%
% INPUT
% control_file control file in MMMx format, must be provided
%
% OUTPUT
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
% entity       output entity at the end of the pipeline, if any, otherwise
%              empty
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
exceptions = {[]};
warnings = 0;
entity = [];

[controller,parse_exceptions] = parse_control_file(control_file);

if ~isempty(parse_exceptions) && ~isempty(parse_exceptions{1})
    exceptions = parse_exceptions;
    return
end

logfid = 1;
close_log = false;
failed = false;
module_exceptions ={};

for module = 1:length(controller)
    fprintf(logfid,'> Executing %s <\n\n',controller(module).name);
    switch lower(controller(module).name)
        case 'logfile'
            close_log = true;
            logfid = fopen(controller(module).options{1},'wt');
        case 'getpdb'
            [entity,module_exceptions] = get_pdb(controller(module).options{1});
            failed = isempty(entity);  
        case 'flex'
            disp('Aber hallo!');
        case 'ensemblefit'
            [entity,module_exceptions,failed] = module_ensemble_fit(controller(module),logfid);
        otherwise
            fprintf(logfid,'Warning: Unknown module %s was ignored.\n\n',controller(module).name);
    end
    for exci = 1:length(module_exceptions)
        if ~isempty(module_exceptions{exci})
            warnings = warnings + 1;
            exceptions{warnings} = module_exceptions{exci};
        end
        if failed
            return
        end
    end
end

if close_log
    fclose(logfid);
end


