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
            failed = isempty(entity); % commands that are essential for further pipeline processing should report failure
        case 'locate'
            module_exceptions = module_locate(controller(module),logfid,entity);
        case 'rigi'
            if ~exist('entity','var')
                fprintf(logfid,'Error: Rigi requires that an entity was loaded that defines rigid bodies.\n\n');
            else
                [entity,exceptions,failed] = module_rigi(controller(module),logfid,entity);
            end
        case 'flex'
            if ~exist('entity','var')
                entity = [];
            end
            [entity,module_exceptions,failed] = module_flex(controller(module),logfid,entity);
        case 'flexrna'
            if ~exist('entity','var')
                entity = [];
            end
            [entity,module_exceptions,failed] = module_flexRNA(controller(module),logfid,entity);
        case 'ensemblefit'
            [entity,module_exceptions,failed] = module_ensemble_fit(controller(module),logfid,entity);
        case 'yasara'
            [entity,module_exceptions,failed] = module_yasara(controller(module),logfid,entity);
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


