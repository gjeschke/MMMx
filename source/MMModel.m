function [exceptions,entity] = MMModel(control_file)
%
% MMModel Run modelling pipeline
%
%   [exceptions,entity] = MMModel(control_file)
%   Parses control file for modelling directives and arguments and runs the
%   specified modelling pipeline
%
% INPUT
% control_file  control file in MMMx format, if missing, an 'open file'
%               window is opened for control file selection
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

my_path = pwd;

if ~exist('control_file','var') || isempty(control_file)
    [control_file,path] = uigetfile('*.mcx');
    cd(path);
end

[controller,parse_exceptions] = parse_control_file(control_file);

[pname,bname,~] = fileparts(control_file);
ext = '.log';
def_logfile = fullfile(pname,[bname ext]);

if ~isempty(parse_exceptions) && ~isempty(parse_exceptions{1})
    exceptions = parse_exceptions;
    cd(my_path);
    return
end

logfid = 1;
close_log = false;
show_report = false;
logfile_name = def_logfile;
failed = false;
module_exceptions ={};
fprintf(logfid,'> Executing %s <\n\n',control_file);

for module = 1:length(controller)
    switch lower(controller(module).name)
        case 'logfile'
            close_log = true;
            if ~isempty(controller(module).options) && ~isempty(controller(module).options{1})
                logfid = fopen(controller(module).options{1},'wt');
                logfile_name = controller(module).options{1};
            else
                logfid = fopen(def_logfile,'wt');
            end
            fprintf(logfid,'> Executing %s <\n\n',control_file);
            timing = tic;
        case 'log'
            logfid = fopen(def_logfile,'wt');
            fprintf(logfid,'> Executing %s <\n\n',control_file);
            timing = tic;
        case 'report'
            show_report = true;
        case 'getpdb'
            [entity,module_exceptions] = get_pdb(controller(module).options{1});
            failed = isempty(entity); % commands that are essential for further pipeline processing should report failure
        case 'locate'
            module_exceptions = module_locate(controller(module),logfid,entity);
        case 'enm'
            [entity,exceptions,failed] = module_enm(controller(module),logfid,entity);
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
        case {'flexrna','flex_rna'}
            if ~exist('entity','var')
                entity = [];
            end
            [entity,module_exceptions,failed] = module_flexRNA(controller(module),logfid,entity);
        case {'ensemblefit','ensemble_fit'}
            [entity,module_exceptions,failed] = module_ensemble_fit(controller(module),logfid,entity);
        case {'yasara','yasararefine','yasara_refine'}
            [entity,module_exceptions,failed] = module_yasara(controller(module),logfid,entity);
        case 'prepare'
            [entity,module_exceptions,failed] = module_prepare(controller(module),logfid);
        case {'ensembleanalysis','ensemble_analysis'}
            [entity,module_exceptions,failed] = module_ensemble_analysis(controller(module),logfid);
        case {'experimentdesign','experiment_design'}
            [entity,module_exceptions,failed] = module_experiment_design(controller(module),logfid);
        case 'visualize'
            module_visualize(controller(module),logfid);
        otherwise
            fprintf(logfid,'Warning: Unknown module %s was ignored.\n\n',controller(module).name);
    end
    for exci = 1:length(module_exceptions)
        if ~isempty(module_exceptions{exci})
            warnings = warnings + 1;
            exceptions{warnings} = module_exceptions{exci};
        end
        if failed
            cd(my_path);
            return
        end
    end
end

runtime = toc(timing);
fprintf(logfid,'\nExecution took ');
hours = floor(runtime/3600);
if hours > 0
    fprintf(logfid,'%i h, ',hours);
end
minutes = floor((runtime-3600*hours)/60);
if minutes > 0
    fprintf(logfid,'%i min, ',minutes);
end
seconds = runtime-3600*hours-60*minutes;
fprintf(logfid,'%4.1f s\n\n',seconds);

fprintf(logfid,'\n*** MMMx modeling finished ***\n');

if close_log
    fclose(logfid);
end

if show_report
    cmd = strtrim(sprintf('notepad %s &\n',logfile_name));
    system(cmd);
end

cd(my_path);

