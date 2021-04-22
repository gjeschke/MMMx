function [entity,exceptions,failed] = module_rigi(control,logfid,entity)
%
% MODULE_RIGI Runs Rigi module
%
%   [entity,exceptions,entity] = MODULE_RIGI(control,logfid,entity)
%   Parses control file for modelling directives and arguments and runs the
%   specified modelling pipeline
%
% INPUT
% control       control structure with fields
%               .name           'rigi', unused
%               .options        options struct
%               .directives     array of struct with directives
%               .entity_output  number of the entity to be output
% logfid        file handle for log file, defaults to console output
% entity        input entity, defines rigid-body coordinates, must exist
%
% OUTPUT
% entity        output entity
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
% failed        flag indicating whether the module failed
%
% SUPPORTED DIRECTIVES
%
% rigid         rigid-body and reference point definition, at least two
%               rigid sections must exist, all other directives can be
%               missing
% separate      can have arguments 'on' or 'off' (default 'on'), separates
%               rigid bodies in entity if 'on'
% resolution    resolution of exhaustive sampling of core restraint space,
%               defaults to 3 Angstroem
% probability   probability threshold, defaults to 0.5
% models        maximum number of models, defaults to 20'000
% maxtime       maximum time (hours), defaults to 48 h
% maxtrials     maximum number of trials, defaults to no limit
% ddr           distance distribution restraints
% plink         peptide link restraints
% rlink         RNA link restraints
% xlink         cross-link restraints
% maxsize       maximum size of the RBA, defaults to 180 Angstroem
% xl_percentage percentage of cross-links that must be fulfilled, defaults
%               to 30%
% parallel      number of trials to be performed in one parallel block,
%               defaults to 10'000
% superimpose   superposition at one of the rigid bodies
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
exceptions = {[]};
warnings = 0;
failed = false;


% set defaults
r_CACA = 3.8; % length between Calpha atoms

save_name = 'MMMx_rigi'; % default name for saving output entity

superimpose = 0;

separate = true;
restraints.resolution = 3;
restraints.p_model = 0.5;
restraints.max_size = 180;
restraints.xlff = 0.3;

if isempty(logfid)
    logfid = 1;
end
options.logfid = logfid;
options.granule = 10000;
options.models = 20000;
options.max_time = 48;


restraints.ddr(1).labels{1} = '';
restraints.xlinks(1).maxr = [];
restraints.xlinks(1).chain1 = '';
restraints.xlinks(1).residue1 = [];
restraints.xlinks(1).chain2 = '';
restraints.xlinks(1).residue2 = [];
restraints.links(1).maxr = [];
restraints.links(1).di = [];
restraints.links(1).type = '';
restraints.links(1).chain1 = '';
restraints.links(1).residue1 = [];
restraints.links(1).chain2 = '';
restraints.links(1).residue2 = [];
restraints.r_bodies(1).chains = [];
restraints.r_bodies(1).ref = [];
restraints.r_bodies(1).indices = [];
restraints.r_bodies(1).points = [];
restraints.r_bodies(1).label = {''};

ddr_poi = 0;
rb_poi = 0;
link_poi = 0;
xlink_poi = 0;
% read restraints
for d = 1:length(control.directives)
    switch lower(control.directives(d).name)
        case 'separate'
            if ~isempty(control.directives(1).options) && strcmpi(control.directives(d).options{1},'off')
                separate = false;
            else
                separate = true;
            end
        case 'save'
            save_name = control.directives(d).options{1};
        case 'maxsize'
            restraints.max_size = str2double(control.directives(d).options{1});
        case 'probability'
            restraints.p_model = str2double(control.directives(d).options{1});
        case 'resolution'
            restraints.resolution = str2double(control.directives(d).options{1});
        case 'xl_percentage'
            restraints.xlff = str2double(control.directives(d).options{1})/100;
        case 'parallel'
            options.granule = str2double(control.directives(d).options{1});
        case 'superimpose'
            superimpose = str2double(control.directives(d).options{1});
        case 'models'
            options.models = round(str2double(control.directives(d).options{1}));
        case 'maxtime'
            options.max_time = str2double(control.directives(d).options{1});
        case 'maxtrials'
            options.max_trials = round(str2double(control.directives(d).options{1}));
        case 'rigid'
            rb_poi = rb_poi + 1;
            chains = '';
            for kc = 1:length(control.directives(d).options)
                chain = control.directives(d).options{kc};
                % remove parentheses if any
                if chain(1) == '('
                    chain = chain(2:end);
                end
                if chain(end) == ')'
                    chain = chain(1:end-1);
                end
                chains = strcat(chains,chain);
                restraints.r_bodies(rb_poi).chains  = chains;
            end
            restraints.r_bodies(rb_poi).ref = zeros(3,5);
            restraints.r_bodies(rb_poi).chain_assign = pad('',3);
            restraints.r_bodies(rb_poi).res_assign = zeros(1,3);
            restraints.r_bodies(rb_poi).labels = cell(1,3);
            % extract reference points
            for kp = 1:3
                adr = control.directives(d).block{kp,1};
                poi1 = strfind(adr,'(');
                poi2 = strfind(adr,')');
                restraints.r_bodies(rb_poi).chain_assign(kp) = adr(poi1+1:poi2-1);
                restraints.r_bodies(rb_poi).res_assign(kp) = str2double(adr(poi2+1:end));
                restraints.r_bodies(rb_poi).labels{kp} = control.directives(d).block{kp,2};
            end
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
                exceptions{warnings} = MException('module_rigi:empty_ddr_block', 'ddr block %i is empty',d);
                record_exception(exceptions{warnings},logfid);
            end
            if args < 3 % misformed restraint line
                warnings = warnings + 1;
                exceptions{warnings} = MException('module_rigi:misformed_ddr', 'ddr restraint has less than three arguments');
                record_exception(exceptions{warnings},logfid);
                failed = true;
                return
            end
            for kr = 1:nr
                restraints.ddr(ddr_poi).site1{kr} = control.directives(d).block{kr,1};
                restraints.ddr(ddr_poi).site2{kr} = control.directives(d).block{kr,2};
                restraints.ddr(ddr_poi).file{kr} = '';
                arg3 = control.directives(d).block{kr,3};
                if arg3(1) == '@'
                    restraints.ddr(ddr_poi).file(kr) = arg3(2:end);
                    restraints.ddr(ddr_poi).r(kr) = [];
                    restraints.ddr(ddr_poi).sigr(kr) = [];
                else
                    restraints.ddr(ddr_poi).r(kr) = str2double(arg3);
                    restraints.ddr(ddr_poi).sigr(kr) = str2double(control.directives(d).block{kr,4});
                end
                if args > 4
                    arg5 = control.directives(d).block{kr,5};
                    if ~isempty(arg5) && arg5(1) == '@'
                        restraints.ddr(ddr_poi).file{kr} = arg5(2:end);
                    end
                end
                for karg = 1:args
                    fprintf(logfid,'  %s',control.directives(d).block{kr,karg});
                end
                fprintf(logfid,'\n');
            end
            fprintf(logfid,'\n\n');
        case 'nlink'
            [nx,~] = size(control.directives(d).block);
            % extract chains, residues, maximum distances
            for kx = 1:nx
                link_poi = link_poi + 1;
                restraints.links(link_poi).type = 'n';
                adr1 = control.directives(d).block{kx,1};
                poi1 = strfind(adr1,'(');
                poi2 = strfind(adr1,')');
                restraints.links(link_poi).chain1 = adr1(poi1+1:poi2-1);
                restraints.links(link_poi).residue1 = str2double(adr1(poi2+1:end));
                adr2 = control.directives(d).block{kx,2};
                poi1 = strfind(adr2,'(');
                poi2 = strfind(adr2,')');
                restraints.links(link_poi).chain2 = adr2(poi1+1:poi2-1);
                restraints.links(link_poi).residue2 = str2double(adr2(poi2+1:end));
                restraints.links(link_poi).di = str2double(control.directives(d).block{kx,3});
                restraints.links(link_poi).maxr = str2double(control.directives(d).block{kx,4});
            end
        case 'plink'
            [nx,~] = size(control.directives(d).block);
            % extract chains, residues, maximum distances
            for kx = 1:nx
                link_poi = link_poi + 1;
                restraints.links(link_poi).type = 'p';
                adr1 = control.directives(d).block{kx,1};
                poi1 = strfind(adr1,'(');
                poi2 = strfind(adr1,')');
                restraints.links(link_poi).chain1 = adr1(poi1+1:poi2-1);
                restraints.links(link_poi).residue1 = str2double(adr1(poi2+1:end));
                adr2 = control.directives(d).block{kx,2};
                poi1 = strfind(adr2,'(');
                poi2 = strfind(adr2,')');
                restraints.links(link_poi).chain2 = adr2(poi1+1:poi2-1);
                restraints.links(link_poi).residue2 = str2double(adr2(poi2+1:end));
                restraints.links(link_poi).di = str2double(control.directives(d).block{kx,3});
                restraints.links(link_poi).maxr = r_CACA*restraints.links(link_poi).di;
            end
        case 'xlink'
            [nx,~] = size(control.directives(d).block);
            % extract chains, residues, maximum distances
            for kx = 1:nx
                xlink_poi = xlink_poi + 1;
                adr1 = control.directives(d).block{kx,1};
                poi1 = strfind(adr1,'(');
                poi2 = strfind(adr1,')');
                restraints.xlinks(xlink_poi).chain1 = adr1(poi1+1:poi2-1);
                restraints.xlinks(xlink_poi).residue1 = str2double(adr1(poi2+1:end));
                adr2 = control.directives(d).block{kx,2};
                poi1 = strfind(adr2,'(');
                poi2 = strfind(adr2,')');
                restraints.xlinks(xlink_poi).chain2 = adr2(poi1+1:poi2-1);
                restraints.xlinks(xlink_poi).residue2 = str2double(adr2(poi2+1:end));
                restraints.xlinks(xlink_poi).maxr = str2double(control.directives(d).block{kx,3});
            end
        otherwise
            warnings = warnings + 1;
            exceptions{warnings} = MException('module_flex:unknown_directive',...
                'directive %s is unknown',lower(control.directives(d).name));
            record_exception(exceptions{warnings},logfid);

    end
end

restraints.xlinks = restraints.xlinks(1:xlink_poi);
restraints.links = restraints.links(1:link_poi);


% separate the rigid bodies if requested, this is the default
if separate
    entity = separate_rigid_bodies(entity,restraints.r_bodies);
end

% Assemble the points defined in each rigid body and make separated
% restraint lists

restraints.points = cell(1,rb_poi);
restraints.pindices = cell(1,rb_poi);
restraints.plabels = cell(1,rb_poi);
restraints.dmat0 = zeros(3*rb_poi);
restraints.auxiliary = zeros(1000,6);

aux_poi = 0;
core_poi = 0;

% Assemble reference points and make intra-rigid-body part of the distance
% matrix
rb_point_assignment = cell(1,rb_poi);
for kr = 1:rb_poi % 
    coor = zeros(500,3);
    point_assignment.chain = pad('',500);
    point_assignment.residue = zeros(1,500);
    point_assignment.points = 3;
    point_assignment.type = cell(1,500);
    for kp = 1:3
        point_assignment.chain(kp) = restraints.r_bodies(kr).chain_assign(kp);
        point_assignment.residue(kp) = restraints.r_bodies(kr).res_assign(kp);
        point_assignment.type{kp} = restraints.r_bodies(kr).labels{kp}; % label
        address = sprintf('(%s)%i',restraints.r_bodies(kr).chain_assign(kp),restraints.r_bodies(kr).res_assign(kp));
        [argsout,entity] = get_label(entity,restraints.r_bodies(kr).labels{kp},'positions',address);
        positions = argsout{1};
        [argsout,entity] = get_label(entity,restraints.r_bodies(kr).labels{kp},'populations',address);
        populations = argsout{1};
        populations = populations/sum(populations);
        coor(kp,:) = populations'*positions;
        restraints.r_bodies(kr).ref(kp,1) = kr;
        restraints.r_bodies(kr).ref(kp,2:4) = coor(kp,:);
    end
    restraints.points{kr} = coor;
    rb_point_assignment{kr} = point_assignment;
    bas = (kr-1)*3;
    for k1 = 1:2
        for k2 = k1:3
            r = norm(coor(k1,:)-coor(k2,:));
            restraints.dmat0(bas+k1,bas+k2) = r;
            restraints.dmat0(bas+k2,bas+k1) = r;
        end
    end
end
restraints.lb = restraints.dmat0; % initialize lower bounds, reference point distances are strictly preserved
restraints.ub = restraints.dmat0; % initialize upper bounds, reference point distances are strictly preserved

% Separate distance restraints into those between reference points,
% auxiliary ones, and delayed restraints that require stemloop attachment,
% add the former ones to the distance matrix

for kblock = 1:ddr_poi
    label1 = restraints.ddr(kblock).labels{1};
    label2 = restraints.ddr(kblock).labels{2};
    for kr = 1:length(restraints.ddr(kblock).site1)
        site1 = restraints.ddr(kblock).site1{kr};
        [chain1,residue1] = split_address(site1);
        [rb1,rp1,pa1] = assign_rb(chain1,residue1,label1,restraints.r_bodies,rb_point_assignment);
        site2 = restraints.ddr(kblock).site2{kr};
        [chain2,residue2] = split_address(site2);
        [rb2,rp2,pa2] = assign_rb(chain2,residue2,label2,restraints.r_bodies,rb_point_assignment);
        if isempty(rb1)
            fprintf(logfid,'Warning: Label %s is not in a rigid body. Restraint will be ignored.\n',site1);
            continue
        end
        if isempty(rb2)
            fprintf(logfid,'Warning: Label %s is not in a rigid body. Restraint will be ignored.\n',site2);
            continue
        end
        if rb1 == rb2
            fprintf(logfid,'Warning: Labels %s and %s are in the same rigid body. Restraint will be ignored.\n',site1,site2);
            continue
        end
        % check whether this is a reference point pair
        if ~isempty(rp1) && ~isempty(rp2)
            fprintf(logfid,'Reference point pair (%i,%i)/(%i,%i) is restrained at %4.1f Angstroem.\n',...
                rb1,rp1,rb2,rp2,restraints.ddr(kblock).r(kr));
            k1 = 3*(rb1-1)+rp1;
            k2 = 3*(rb2-1)+rp2;
            restraints.dmat0(k1,k2) = restraints.ddr(kblock).r(kr);
            restraints.dmat0(k2,k1) = restraints.ddr(kblock).r(kr);
            restraints.lb(k1,k2) = restraints.ddr(kblock).r(kr) - restraints.ddr(kblock).sigr(kr);
            restraints.lb(k2,k1) = restraints.ddr(kblock).r(kr) - restraints.ddr(kblock).sigr(kr);
            restraints.ub(k1,k2) = restraints.ddr(kblock).r(kr) + restraints.ddr(kblock).sigr(kr);
            restraints.ub(k2,k1) = restraints.ddr(kblock).r(kr) + restraints.ddr(kblock).sigr(kr);
            core_poi = core_poi + 1;
            restraints.core(core_poi,:) = [rb1 rp1 rb2 rp2 restraints.ddr(kblock).r(kr) restraints.ddr(kblock).sigr(kr)];
        else % this is an auxiliary restraint
            % check whether new points need to be recorded
            if isempty(pa1)
                point_assignment = rb_point_assignment{rb1};
                points = point_assignment.points + 1;
                pa1 = points;
                point_assignment.points = points;
                point_assignment.chain(points) = chain1;
                point_assignment.residue(points) = residue1;
                point_assignment.type{points} = label1; % label
                rb_point_assignment{rb1} = point_assignment;
                coor = restraints.points{rb1};
                address = sprintf('(%s)%i',chain1,residue1);
                [argsout,entity] = get_label(entity,label1,'positions',address);
                positions = argsout{1};
                [argsout,entity] = get_label(entity,label1,'populations',address);
                populations = argsout{1};
                populations = populations/sum(populations);
                coor(points,:) = populations'*positions;
                restraints.points{rb1} = coor;
            end
            if isempty(pa2)
                point_assignment = rb_point_assignment{rb2};
                points = point_assignment.points + 1;
                pa2 = points;
                point_assignment.points = points;
                point_assignment.chain(points) = chain2;
                point_assignment.residue(points) = residue2;
                point_assignment.type{points} = label2; % label
                rb_point_assignment{rb2} = point_assignment;
                coor = restraints.points{rb2};
                address = sprintf('(%s)%i',chain2,residue2);
                [argsout,entity] = get_label(entity,label2,'positions',address);
                positions = argsout{1};
                [argsout,entity] = get_label(entity,label2,'populations',address);
                populations = argsout{1};
                populations = populations/sum(populations);
                coor(points,:) = populations'*positions;
                restraints.points{rb2} = coor;
            end
            aux_poi = aux_poi + 1;
            restraints.auxiliary(aux_poi,:) = [rb1 pa1 rb2 pa2 restraints.ddr(kblock).r(kr) restraints.ddr(kblock).sigr(kr)];
        end
    end
end
restraints.auxiliary = restraints.auxiliary(1:aux_poi,:);

% determine points for crosslinks
for kx = 1:length(restraints.xlinks)
    [rb1,~,pa1] = assign_rb(restraints.xlinks(kx).chain1,restraints.xlinks(kx).residue1,...
        'CA',restraints.r_bodies,rb_point_assignment);
    [rb2,~,pa2] = assign_rb(restraints.xlinks(kx).chain2,restraints.xlinks(kx).residue2,...
        'CA',restraints.r_bodies,rb_point_assignment);
    if isempty(pa1)
        point_assignment = rb_point_assignment{rb1};
        points = point_assignment.points + 1;
        pa1 = points;
        point_assignment.points = points;
        point_assignment.chain(points) = restraints.xlinks(kx).chain1;
        point_assignment.residue(points) = restraints.xlinks(kx).residue1;
        point_assignment.type{points} = 'CA'; 
        rb_point_assignment{rb1} = point_assignment;
        resstr = sprintf('R%i',restraints.xlinks(kx).residue1);
        index = entity.(restraints.xlinks(kx).chain1).(resstr).CA.tab_indices;
        coor = restraints.points{rb1};
        coor(points,:) = entity.xyz(index,:);
        restraints.points{rb1} = coor;
    end
    if isempty(pa2)
        point_assignment = rb_point_assignment{rb2};
        points = point_assignment.points + 1;
        pa2 = points;
        point_assignment.points = points;
        point_assignment.chain(points) = restraints.xlinks(kx).chain2;
        point_assignment.residue(points) = restraints.xlinks(kx).residue2;
        rb_point_assignment{rb2} = point_assignment;
        coor = restraints.points{rb2};
        resstr = sprintf('R%i',restraints.xlinks(kx).residue2);
        index = entity.(restraints.xlinks(kx).chain2).(resstr).CA.tab_indices;
        coor(points,:) = entity.xyz(index,:);
        restraints.points{rb2} = coor;
    end
    restraints.xlinks(kx).ref_indices = [rb1,pa1,rb2,pa2];
end

% determine anchor points for flexible linkers
for kl = 1:length(restraints.links)
    switch restraints.links(kl).type
        case 'p'
            site = 'CA';
        case 'n'
            site = 'C5_';
    end      
    [rb1,~,pa1] = assign_rb(restraints.links(kl).chain1,restraints.links(kl).residue1,...
        site,restraints.r_bodies,rb_point_assignment);
    [rb2,~,pa2] = assign_rb(restraints.links(kl).chain2,restraints.links(kl).residue2,...
        site,restraints.r_bodies,rb_point_assignment);
    if isempty(pa1)
        point_assignment = rb_point_assignment{rb1};
        points = point_assignment.points + 1;
        pa1 = points;
        point_assignment.points = points;
        point_assignment.chain(points) = restraints.links(kl).chain1;
        point_assignment.residue(points) = restraints.links(kl).residue1;
        point_assignment.type{points} = site;
        rb_point_assignment{rb1} = point_assignment;
        resstr = sprintf('R%i',restraints.links(kl).residue1);
        index = entity.(restraints.links(kl).chain1).(resstr).(site).tab_indices;
        coor = restraints.points{rb1};
        coor(points,:) = entity.xyz(index,:);
        restraints.points{rb1} = coor;
    end
    if isempty(pa2)
        point_assignment = rb_point_assignment{rb2};
        points = point_assignment.points + 1;
        pa2 = points;
        point_assignment.points = points;
        point_assignment.chain(points) = restraints.links(kl).chain2;
        point_assignment.residue(points) = restraints.links(kl).residue2;
        point_assignment.type{points} = site;
        rb_point_assignment{rb2} = point_assignment;
        coor = restraints.points{rb2};
        resstr = sprintf('R%i',restraints.links(kl).residue2);
        index = entity.(restraints.links(kl).chain2).(resstr).(site).tab_indices;
        coor(points,:) = entity.xyz(index,:);
        restraints.points{rb2} = coor;
    end
    restraints.links(kl).ref_indices = [rb1,pa1,rb2,pa2];
end

% restrict point arrays for rigid bodies to actual numbers of points
for kb = 1:rb_poi
    points = restraints.points{kb};
    point_assignment = rb_point_assignment{kb};
    restraints.points{kb} = points(1:point_assignment.points,:);
    point_assignment.chain = point_assignment.chain(1:point_assignment.points);
    point_assignment.residue = point_assignment.residue(1:point_assignment.points);
    point_assignment.type = point_assignment.type(1:point_assignment.points);
    rb_point_assignment{kb} = point_assignment;
end
% store point assignment, although it is not used in rigi_engine
restraints.rb_point_assignment = rb_point_assignment;

[transforms,diagnostics] = rigi_engine(entity,restraints,options);

% superimpose models by correcting translation/rotation of rigid bodies
if superimpose > 0
    if superimpose > rb_poi
        fprintf(logfid,'Superposition on rigid body %i failed as there are only %i rigid bodies.\n',superimpose,n_rb); 
    else
        [n_rba,~] = size(transforms);
        bas = 6*(superimpose-1);
        for krba = 1:n_rba
            % separate translation and rotation matrices for reference
            % rigid body
            transrot0 = transforms(krba,bas+1:bas+6);
            rotmat = transrot2affine(zeros(1,3),transrot0(4:6));
            % invert rotation
            rotmat(1:3,1:3) = rotmat(1:3,1:3)';
            % the following is the inverse translation
            transmat = transrot2affine(-transrot0(1:3),zeros(1,3));
            for rb = 1:rb_poi
                bas2 = 6*(rb-1);
                transrot = transforms(krba,bas2+1:bas2+6);
                % original affine transformation matrix
                origmat = transrot2affine(transrot(1:3),transrot(4:6));
                newmat = rotmat*transmat*origmat;
                [trans,euler] = affine2transrot(newmat);
                transforms(krba,bas2+1:bas2+3) = trans;
                transforms(krba,bas2+4:bas2+6) = euler;
            end
        end
    end
end

% write diagnostics information into logfile

hours = floor(diagnostics.runtime/3600);
minutes = floor((diagnostics.runtime-3600*hours)/60);
seconds = round(diagnostics.runtime - 3600*hours - 60*minutes);

done = 'not';
if diagnostics.completed
    done = '';
end

fprintf(logfid,'\nRigi run was %s completed in %i h %i min %i s\n\n',done,hours,minutes,seconds);
fprintf(logfid,'%i models were generated in %i trials\n',diagnostics.success,diagnostics.trials);
if diagnostics.solutions > diagnostics.success 
    fprintf(logfid,'They were selected from %i initials solutions by clustering with resolution %4.1f Angstroem.\n',...
        diagnostics.solutions,diagnostics.cluster_resolution);
    hours = floor(diagnostics.cluster_time/3600);
    minutes = floor((diagnostics.cluster_time-3600*hours)/60);
    seconds = round(diagnostics.cluster_time - 3600*hours - 60*minutes);
    fprintf(logfid,'Clustering took %i h %i min %i s\n',hours,minutes,seconds);
end
fprintf(logfid,'The worst encountered resolution was %4.1f Angstroem\n\n',diagnostics.resolution);

% store entity in MMMx|RigiFlex format
entity.rba = transforms;
entity.rba_populations = diagnostics.probabilities/sum(diagnostics.probabilities);

% add information on rigid-body assignment
entity.rigidbodies = cell(1,length(restraints.r_bodies));
entity.rba_chains = cell(1,length(restraints.r_bodies));
[m,~] = size(entity.xyz);
for kr = 1:length(restraints.r_bodies)
    chains = restraints.r_bodies(kr).chains;
    entity.rba_chains{kr} = chains;
    all_indices = zeros(m,1);
    atnum = 0;
    for kc = 1:length(chains)
        entity.(chains(kc)).rigidbody = kr;
        address = sprintf('(%s)',restraints.r_bodies(kr).chains(kc));
        [~,indices] = get_coor(entity,address);
        mc = length(indices);
        all_indices(1+atnum:mc+atnum) = indices;
        atnum = atnum + mc;
    end
    entity.rigidbodies{kr} = all_indices(1:atnum);
end

% save generated entity
save(save_name,'entity');

function record_exception(exception,logfid)

fprintf(logfid,'### rigi exception %s ###\n',exception.message);

function [chain,residue] = split_address(address)

poi1 = strfind(address,'(');
poi2 = strfind(address,')');
chain = address(poi1+1:poi2-1);
residue = str2double(address(poi2+1:end));

function [rb,rp,pa] = assign_rb(chain,residue,type,rigid_bodies,rb_point_assignment)
% assign a site to a rigid body and reference point within a rigid body

rb = [];
rp = [];
pa = [];
for kb = 1:length(rigid_bodies)
    % reference points take precedence
    if contains(rigid_bodies(kb).chain_assign,chain)
        rb = kb;
        for kp = 1:3
            if chain == rigid_bodies(kb).chain_assign(kp) && ...
                    residue == rigid_bodies(kb).res_assign(kp)
                rp = kp;
            end
        end
        point_assignment = rb_point_assignment{kb};
        for kp = 1:point_assignment.points
            if chain == point_assignment.chain(kp) && ...
                    residue == point_assignment.residue(kp) && ...
                    strcmpi(type,point_assignment.type{kp})
                pa = kp;
            end
        end
    elseif contains(rigid_bodies(kb).chains,chain)
        rb = kb;
        point_assignment = rb_point_assignment{kb};
        for kp = 1:point_assignment.points
            if chain == point_assignment.chain(kp) && ...
                    residue == point_assignment.residue(kp) && ...
                    strcmpi(type,point_assignment.type{kp})
                pa = kp;
            end
        end
    end
end