function [exceptions,failed] = module_locate(control,logfid,entity)
%
% MODULE_LOCATE Runs module for locating a site by distance restraints
%
%   [entity,exceptions,entity] = MODULE_LOCATE(control,logfid,entity)
%
% INPUT
% control       control structure with fields
%               .name           'locate', unused
%               .options        options struct
%               .directives     array of struct with directives
%               .entity_output  number of the entity to be output
% logfid        file handle for log file, defaults to console output
% entity        input entity, defines rigid-body coordinates, must exist
%
% OUTPUT
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
% failed        flag indicating whether the module failed
%
% SUPPORTED DIRECTIVES
%
% reference     reference (GPS-like) distance restraints between sites in
%               the template entity and the site to be localized
% save          name for the cube file of the output isodensity surface, if
%               not specified, name is MMMx_locate.mrc
% cube          side length of density cube, defaults to 75 Angstroem
% grid          number of grid points along one side of the density cube, 
%               defaults to 176
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% initialize empty output
exceptions = {[]};
warnings = 0;
failed = false;


% set defaults

save_name = 'MMMx_locate.mrc'; % default name for saving output entity

restraints.p_model = 0.5;
if ~isempty(control.options)
    restraints.p_model = str2double(control.options{1});
end

if isempty(logfid)
    logfid = 1;
end


restraints.reference(1).label = '';

ref_poi = 0;
% read restraints
for d = 1:length(control.directives)
    switch lower(control.directives(d).name)
        case 'getpdb'
            entity = get_pdb(control.directives(d).options{1});
        case 'getalphafold'
            entity = get_AF(control.directives(d).options{1});
        case 'reference'
            ref_poi = ref_poi + 1; % increase reference block counter
            restraints.reference(ref_poi).label = control.directives(d).options{1};
            [nr,args] = size(control.directives(d).block);
            if nr > 0 % at least one restraint in this block
                restraints.reference(ref_poi).site{nr} = ' ';
                restraints.reference(ref_poi).r(nr) = 0;
                restraints.reference(ref_poi).sigr(nr) = 0;
                restraints.reference(ref_poi).file{nr} = '*';
            else % empty block
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_locate:empty_reference_block', 'reference block %i is empty',d);
                record_exception(exceptions{warnings},logfid);
            end
            if args < 2 % misformed restraint line
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_locate:misformed_reference_restraint', 'reference restraints have less than two arguments');
                record_exception(exceptions{warnings},logfid);
                failed = true;
                return
            end
            for kr = 1:nr
                restraints.reference(ref_poi).site{kr} = control.directives(d).block{kr,1};
                restraints.reference(ref_poi).file{kr} = '';
                arg2 = control.directives(d).block{kr,2};
                if arg2(1) == '@'
                    restraints.reference(ref_poi).file{kr} = arg2(2:end);
                    restraints.reference(ref_poi).r(kr) = NaN;
                    restraints.reference(ref_poi).sigr(kr) = NaN;
                else
                    restraints.reference(ref_poi).r(kr) = str2double(arg2);
                    restraints.reference(ref_poi).sigr(kr) = str2double(control.directives(d).block{kr,3});
                end
                if args > 3
                    arg4 = control.directives(d).block{kr,4};
                    if ~isempty(arg4) && arg4(1) == '@'
                        restraints.reference(ref_poi).file{kr} = arg4(2:end);
                    end
                end
            end
        case 'save'
            save_name = control.directives(d).options{1};
            % add extension if not present
            if ~contains(save_name,'.')
                save_name = strcat(save_name,'.mrc');
            end
        case 'cube'
            restraints.cube = str2double(control.directives(d).options{1});
        case 'grid'
            restraints.grid = round(str2double(control.directives(d).options{1}));
        otherwise
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_locate:unknown_directive',...
                'directive %s is unknown',lower(control.directives(d).name));
            record_exception(exceptions{warnings},logfid);
    end
end

% count reference restraints
ref_num = 0;
for kr = 1:ref_poi
    ref_num = ref_num + length(restraints.reference(kr).site);
end

% location fails if there are less than three reference restraints
if ref_num < 3
    warnings = warnings + 1;
    exceptions{warnings} = MException('module_locate:too_few_reference_points',...
        'ERROR: Locate requires at least 3 reference points, only %i are provided (Aborting)',ref_num);
    record_exception(exceptions{warnings},logfid);
    failed = true;
    return
end
% with only three reference restraints, location is ambiguous
if ref_num == 3
    warnings = warnings + 1;
    exceptions{warnings} = MException('module_locate:ambiguous_location',...
        'Warning: Location is ambiguous with only three reference points');
    record_exception(exceptions{warnings},logfid);
end

% find and store mean spin label coordinates for reference points
for kr = 1:length(restraints.reference)
    for nr = 1:length(restraints.reference(kr).site)
        [argsout,entity] = get_label(entity,restraints.reference(kr).label,'positions',restraints.reference(kr).site{nr});
        if isempty(argsout)
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_locate:reference_point_missing',...
                'ERROR: Reference site %s missing or cannot be labelled (Aborting)',restraints.reference(kr).site{nr});
            record_exception(exceptions{warnings},logfid);
            failed = true;
            return
        end
        positions = argsout{1};
        [argsout,entity] = get_label(entity,restraints.reference(kr).label,'populations',restraints.reference(kr).site{nr});
        if isempty(argsout) || isempty(positions) || isempty(argsout{1})
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_locate:reference_point_cannot_be_labelled',...
                'ERROR: Reference site %s exists, but cannot be labelled (Aborting)',restraints.reference(kr).site{nr});
            record_exception(exceptions{warnings},logfid);
            failed = true;
            return
        end
        populations = argsout{1};
        populations = populations/sum(populations);
        restraints.reference(kr).coor{nr} = populations'*positions;
        if ~isempty(restraints.reference(kr).file{nr}) % if a distribution file exists, load distribution
            exp_data = load_distance_distribution(restraints.reference(kr).file{nr});
            rax = exp_data(:,1)';
            restraints.reference(kr).rax{nr} = rax; % conversion to Angstroem
            distr = exp_data(:,2)';
            distr = distr/sum(distr);
            restraints.reference(kr).distr{nr} = distr;
            % set mean distance and standard deviation, if not provided
            if isnan(restraints.reference(kr).r(nr))
                rmean = sum(rax.*distr);
                restraints.reference(kr).r(nr) = rmean;
                restraints.reference(kr).sigr(nr) = sqrt(sum(distr.*(rax-rmean).^2));
            end
        else % if not provided, make distributions
            rax = 0:0.05:restraints.reference(kr).r(nr)+3*restraints.reference(kr).sigr(nr);
            distr = exp(-(rax-restraints.reference(kr).r(nr)).^2/(2*restraints.reference(kr).sigr(nr)^2));
            distr = distr/sum(distr);
            restraints.reference(kr).rax{nr} = rax;
            restraints.reference(kr).distr{nr} = distr;
        end
    end
end

restraints.save_name = save_name;

rigid_multilateration(restraints,logfid);


function record_exception(exception,logfile)

fprintf(logfile,'### locate exception: %s ###\n',exception.message);