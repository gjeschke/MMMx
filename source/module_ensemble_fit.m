function [entity,exceptions,failed,restraints] = module_ensemble_fit(control,logfid)
%
% MMModel Run modelling pipeline
%
%   [exceptions,entity] = MMModel(control_file)
%   Parses control file for modelling directives and arguments and runs the
%   specified modelling pipeline
%
% INPUT
% control       control structure with fields
%               .name           'ensemblefit', unused
%               .options        options struct
%               .directives     array of struct with directives
%               .entity_output  number of the entity to be output
% logfid        file handle for log file
%
% OUTPUT
% entity        output entity, empty in NOFIT mode
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
% failed        flag indicated whther the module failed
% restraints    restraint structure with fit values
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
exceptions = {[]};
warnings = 0;
entity = [];
failed = false;

restraints.ddr(1).labels{1} = '';
restraints.saxs(1).data{1} = '';
restraints.sans(1).data{1} = '';

ddr_poi = 0;
saxs_poi = 0;
sans_poi = 0;

% read restraints
for d = 1:length(control.directives)
    switch control.directives(d).name
        case 'ddr'
            ddr_poi = ddr_poi + 1; % increase ddr block counter
            fprintf(logfid,'ddr %s',control.directives(d).options{1}); % echo directive to log file
            restraints.ddr(ddr_poi).labels{1} = control.directives(d).options{1};
            if length(control.directives(d).options) > 1 % different labels
                restraints.ddr(ddr_poi).labels{2} = control.directives(d).options{2};
                fprintf(logfid,' %s\n',control.directives(d).options{2});
            else % same label at both sites
                restraints.ddr(ddr_poi).labels{2} = control.directives(d).options{1};
                fprintf(logfid,'\n');
            end
            [nr,args] = size(control.directives(d).block);
            if nr > 0 % at least one restraint in this block
                restraints.ddr(ddr_poi).site1{nr} = ' ';
                restraints.ddr(ddr_poi).site2{nr} = ' ';
                restraints.ddr(ddr_poi).r(nr) = 0;
                restraints.ddr(ddr_poi).sigr(nr) = 0;
                restraints.ddr(ddr_poi).file{nr} = '*';
            else % empty block
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_ensemble_fit:empty_ddr_block', 'ddr block %i is empty',d);
            end
            if args < 3 % misformed restraint line
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_ensemble_fit:misformed_ddr', 'ddr restraint has less than three arguments');
                failed = true;
                return
            end
            for kr = 1:nr
                restraints.ddr(ddr_poi).site1{kr} = control.directives(d).block{kr,1};
                restraints.ddr(ddr_poi).site2{kr} = control.directives(d).block{kr,2};
                arg3 = control.directives(d).block{kr,3};
                if arg3(1) == '@'
                    restraints.ddr(ddr_poi).file(kr) = arg3(2:end);
                else
                    restraints.ddr(ddr_poi).r(kr) = str2double(arg3);
                    restraints.ddr(ddr_poi).sigr(kr) = str2double(control.directives(d).block{kr,4});
                end
                if args > 4
                    arg5 = control.directives(d).block{kr,5};
                    if arg5(1) == '@'
                        restraints.ddr(ddr_poi).file{kr} = arg5(2:end);
                    end
                end
                for karg = 1:args
                    fprintf(logfid,'  %s',control.directives(d).block{kr,karg});
                end
                fprintf(logfid,'\n');
            end
            fprintf(logfid,'\n\n');
    end
end

restraints.ddr = restraints.ddr(1:ddr_poi);
restraints.saxs = restraints.saxs(1:saxs_poi);
restraints.sans = restraints.sans(1:sans_poi);

