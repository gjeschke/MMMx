function module_visualize(control,logfid)
%
% MODULE_VISUALIZE    Visualizes MMMx results by calling ChimeraX
%
%   MODULE_VISUALIZE(control,logfid)
%   Visualization functions that are specific to ensemble modelling and
%   to the use of spin labels
%
% INPUT
% control       control structure with fields
%               .name           'visualize', unused
%               .options        options struct
%               .directives     array of struct with directives
%               .entity_output  number of the entity to be output
% logfid        file handle for log file, defaults to console output
%
% SUPPORTED DIRECTIVES (an ordered list is processed)
%
% addpdb        add conformers by reading pdb files
% getens        get ensemble file
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

global hMain % if MMM is open, these are handles to the MMM gui

% set defaults

figure_format = 'png';
script = 'MMMx.mmm';
execution = false;
normalize = true;

% default output to Matlab console if no log file identifiere was provided 
if isempty(logfid)
    logfid = 1;
end

commands = cell(1,1000); % command list
densities = cell(1,1000); % command list
graphics = cell(1,1000); % graphics list

cmd_poi = 0; % command pointer
fig_poi = 0;
dens_poi = 0;
% reorganize command line arguments
for d = 1:length(control.directives)
    clear cmd
    cmd.name = lower(control.directives(d).name);
    switch lower(control.directives(d).name)
        case {'getens','input'}
            cmd.input = control.directives(d).options{1};
            % allow for input of zipped ensembles
            [~,~,extension] = fileparts(cmd.input);
            if strcmpi(extension,'.zip')
                filenames = unzip(cmd.input);
                for f = 1:length(filenames)
                    [~,~,ext] = fileparts(filenames{f});
                    if strcmpi(ext,'.ens')
                        cmd.input = filenames{f};
                    end
                end
            end
            entity = get_ensemble(cmd.input);
            args = split(control.directives(d).options{1},'.');
            fname = args{1};
            svname = sprintf('MMMx_visualize_%s.pdb',fname);
            put_pdb(entity,svname);
            all_pdb(1).name = svname;
            pop = entity.populations;
        case {'getpdb','import'}
            entity = get_pdb(control.directives(d).options{1});
            fname = sprintf('MMMx_visualize_%s.pdb',control.directives(d).options{1});
            put_pdb(entity,fname);
            all_pdb(1).name = fname;
            pop = 1;
        case {'getalphafold'}
            entity = get_AF(control.directives(d).options{1});
            fname = sprintf('MMMx_visualize_%s.pdb',control.directives(d).options{1});
            put_pdb(entity,fname);
            all_pdb(1).name = fname;
            pop = 1;
        case {'get_zenodo'}
            args = split(control.directives(d).options{1},'.');
            fname = args{2};
            svname = sprintf('MMMx_visualize_%s.pdb',fname);
            k = 2;
            while k < length(args)
                k = k + 1;
                fname = sprintf('%s.%s',fname,args{k});
            end
            entity = get_zenodo(args{1},fname);
            put_pdb(entity,svname);
            all_pdb(1).name = svname;
            pop = entity.populations;
        case 'addpdb'
            all_pdb = dir(control.directives(d).options{1});
            pop = ones(length(all_pdb),1)/length(all_pdb);
        case 'figures'
            figure_format = control.directives(d).options{1};
        case 'normalize'
            if strcmpi(control.directives(d).options{1},'off')
                normalize = false;
            end
        case 'color'
            cmd_poi = cmd_poi + 1;
            cmd.address = control.directives(d).options{1};
            if length(control.directives(d).options) >= 4
                cmd.rgb = [str2double(control.directives(d).options{2}),str2double(control.directives(d).options{3}),str2double(control.directives(d).options{4})];
            else
                cmd.rgb = get_svg_color(control.directives(d).options{2});
            end
            commands{cmd_poi} = cmd;
        case 'colorscheme'
            cmd_poi = cmd_poi + 1;
            cmd.address = control.directives(d).options{1};
            cmd.scheme = '';
            for arg = 2:length(control.directives(d).options)
                cmd.scheme = sprintf('%s %s',cmd.scheme,control.directives(d).options{arg});
            end
            commands{cmd_poi} = cmd;
        case 'label'
            cmd_poi = cmd_poi + 1;
            cmd.address = control.directives(d).options{1};
            cmd.type = 'mtsl';
            if length(control.directives(d).options) >= 2
                cmd.type = control.directives(d).options{2};
            end
            commands{cmd_poi} = cmd;
        case 'show'
            cmd_poi = cmd_poi + 1;
            cmd.address = control.directives(d).options{1};
            cmd.mode = control.directives(d).options{2};
            commands{cmd_poi} = cmd;
        case 'density'
            dens_poi = dens_poi + 1;
            cmd.fname = control.directives(d).options{1};  
            cmd.level = '';
            cmd.opacity = '';
            cmd.rgb = [0.75,0,0];
            if length(control.directives(d).options) > 1
                cmd.level = control.directives(d).options{2};
            end
            if length(control.directives(d).options) > 2
                cmd.opacity = control.directives(d).options{3};
            end
            if length(control.directives(d).options) >= 6
                cmd.rgb = [str2double(control.directives(d).options{4}),str2double(control.directives(d).options{5}),str2double(control.directives(d).options{6})];
            else
                cmd.rgb = get_svg_color(control.directives(d).options{4});
            end
            densities{dens_poi} = cmd;
        case 'graphics'
            fig_poi = fig_poi + 1;
            cmd.fname = '';
            cmd.view = '';
            if ~isempty(control.directives(d).options)
                cmd.fname = control.directives(d).options{1};                
                cmd.mode = figure_format;
                if length(control.directives(d).options) > 1
                    cmd.mode = control.directives(d).options{2};
                    cmd.view = '';
                    for k = 3:length(control.directives(d).options)
                        cmd.view = sprintf('%s %s',cmd.view,control.directives(d).options{k});
                    end
                end
            end
            graphics{fig_poi} = cmd;
        case 'script'
            script = control.directives(d).options{1};
            [pname,fname,ext] = fileparts(script);
            if isempty(ext)
                ext = '.mmm';
            end
            script = fullfile(pname,strcat(fname,ext));
        case 'execute'
            execution = true;
        case 'isosurface'
            cmd_poi = cmd_poi + 1;
            cmd.density = control.directives(d).options{1};
            cmd.options.colorscheme = '';
            if length(control.directives(d).options) > 1 % property file exists
                cmd.property = control.directives(d).options{2};
                cmd.options.colorscheme = 'electrostatic';
            else
                cmd.property = ''; 
            end
            cmd.options.level = 0.999;
            cmd.options.camvec = [1,0,0];
            cmd.options.opaqueness = 1;
            cmd.options.limits = [];
            cmd.options.camupvec = [0,1,0];
            cmd.options.figname = 'isosurface.png';
            [n,~] = size(control.directives(d).block);
            % set requested options
            for k = 1:n
                switch lower(control.directives(d).block{k,1})
                    case 'level'
                        cmd.options.level = str2double(control.directives(d).block{k,2});
                    case 'camvec'
                        cmd.options.camvec(1) = str2double(control.directives(d).block{k,2});
                        cmd.options.camvec(2) = str2double(control.directives(d).block{k,3});
                        cmd.options.camvec(3) = str2double(control.directives(d).block{k,4});
                    case 'camupvec'
                        cmd.options.camupvec(1) = str2double(control.directives(d).block{k,2});
                        cmd.options.camupvec(2) = str2double(control.directives(d).block{k,3});
                        cmd.options.camupvec(3) = str2double(control.directives(d).block{k,4});
                    case 'limits'
                        if strcmpi(cmd.options.limits,'adapted')
                            cmd.options.limits = NaN;
                        else
                            cmd.options.limits = str2double(control.directives(d).block{k,2});
                        end
                    case 'figname'
                        cmd.options.figname = control.directives(d).block{k,2};
                    case 'opaqueness'
                        cmd.options.opaqueness = str2double(control.directives(d).block{k,2});
                    case 'colorscheme'
                        cmd.options.colorscheme = control.directives(d).block{k,2};
                end
            end
            commands{cmd_poi} = cmd;            
        otherwise
            fprintf(logfid,'directive %s is unknown and will be ignored',lower(control.directives(d).name));
    end
end
all_tags = ':';
if exist('all_pdb','var') 
    for k = 1:length(all_pdb)
        fid=fopen(all_pdb(k).name);
        tline = fgetl(fid);
        fclose(fid);
        if length(tline)>=66
            idCode=tline(63:66);
            if ~strcmpi(strtrim(idCode),idCode), idCode = ''; end
        else
            idCode='';
        end
        stag = idCode;
        id=tag2id(stag,all_tags);
        poi=1;
        while ~isempty(id)
            stag = sprintf('%s_%i',idCode,poi);
            poi=poi+1;
            id = tag2id(stag,all_tags);
        end
        all_tags=sprintf('%s%s:',all_tags,stag);
    end
    if normalize
        pop = pop/max(pop);
    end
end

graphics = graphics(1:fig_poi);
densities = densities(1:dens_poi);
commands = commands(1:cmd_poi);


% run the command list for all ensemble members, this creates the
% visualization

if ~exist('all_pdb','var') || isempty(all_pdb)
    for c = 1:cmd_poi
        cmd = commands{c};
        switch cmd.name
            case {'isosurface'}
                visualize_isosurface(cmd.density,cmd.property,cmd.options);
        end
    end
    return
end

ofid = fopen(script,'wt');
fprintf(ofid,'%% MMMx visualization script\n');
fprintf(ofid,'new !\n');
pop_encode_transparency = true;
for k = 1:length(all_pdb)
    tag = id2tag(k,all_tags);
    sadr = sprintf('[%s]{:}',tag);
    fprintf(ofid,'pdbload %s\n',all_pdb(k).name);
    for c = 1:cmd_poi
        cmd = commands{c};
        switch cmd.name
            case {'isosurface'}
                if k == 1
                    visualize_isosurface(cmd.density,cmd.property,cmd.options);
                end
            case {'show'}
                if strcmpi(cmd.mode,'snake')
                    if length(all_pdb) == 1
                        for conformer = 1:length(pop)
                            show_snake(ofid,pop(conformer),sprintf('[%s]{%i}',tag,conformer),cmd.address);
                        end
                    else
                        show_snake(ofid,pop(k),sadr,cmd.address);
                    end
                    pop_encode_transparency = false;
                else
                    fprintf(ofid,'show %s %s\n',strcat(sadr,cmd.address),cmd.mode);
                end
            case {'color'}
                fprintf(ofid,'color %s %6.3f%6.3f%6.3f\n',strcat(sadr,cmd.address),cmd.rgb);
            case {'colorscheme'}
                fprintf(ofid,'colorscheme %s %s\n',strcat(sadr,cmd.address),cmd.scheme);
            case {'label'}
                fprintf(ofid,'rotamers %s %s ambient\n',strcat(sadr,cmd.address),cmd.type);
                fprintf(ofid,'label %s %s ambient\n',strcat(sadr,cmd.address),cmd.type);
        end
    end
end
% transparency slows down drawing. Hence, it is applied only after
% everything else was drawn
if length(all_pdb) > 1
    for k = 1:length(all_pdb)
        tag = id2tag(k,all_tags);
        sadr = sprintf('[%s]',tag);
        if pop_encode_transparency
            fprintf(ofid,'transparency %s(:) %5.3f\n',sadr,pop(k));
        end
    end
else
    tag = id2tag(1,all_tags);
    sadr = sprintf('[%s]',tag);
    for k = 1:length(pop)
        if pop_encode_transparency
            fprintf(ofid,'transparency %s{%i}(:) %5.3f\n',sadr,k,pop(k));
        end
    end
end
% generate all specified density isosurfaces
for d = 1:dens_poi
    cmd = densities{d};
    fprintf(ofid,'density %s %s %s %5.3f %5.3f %5.3f\n',cmd.fname,cmd.level,cmd.opacity,cmd.rgb);
end

% make all specified graphics (files)
for g = 1:fig_poi
    cmd = graphics{g};
    if ~isempty(cmd.view)
        fprintf(ofid,'view %s\n',cmd.view);
    end
    fprintf(ofid,'detach\n');
    fprintf(ofid,'zoom out\n');
    if isempty(cmd.fname) % copy to clipboard
        fprintf(ofid,'copy\n');
    else
        fprintf(ofid,'copy %s %s\n',cmd.fname,cmd.mode);
    end
end

fclose(ofid);

if execution
    script_file = fullfile(pwd,script);
    MMM_prototype('execute_script_callback',hMain.button_script,script_file,hMain);
end

function tag = id2tag(id,tags,codelist,delimiter)
% function tag=id2tag(id,tags,codelist,delimiter)
%
% Returns the string tag corresponding to an identification code from
% a list of possible tags given in a string with colon (:)
% separation, another separator can be selected by optional parameter
% delimiter
% if an array codelist is given the identifier is of the same type as the
% elements in codelist, the sequence in the codelist must correspond to the
% one in the tags string
% if codelist is missing, the id is an integer corresponding to the
% position of the tag in string tags
% if the id is larger than the number of tags in the tags string an empty
% tag is given back
% if the id appears twice or more in the code list, the first tag is given
% back
%
% id        identifier, for example 7
% tags      (string) tag list, example ':H:He:Li:Be:C:N:O:F:Ne:'
% codelist  [optional] (array of ids), example [1,4,7,9,12,14,16,19,20]
% delimiter [optional] (char), alternative delimiter, defaults to ':'
%
% tag       (string) tag
%           if codelist is missing, tag for the example will be 'O'
%           if codelist is present, tag for the example will be 'Li'
%
% G. Jeschke, 2009

tag=''; % empty default output
if nargin < 4
    delimiter=':'; % colon as default delimiter
end

if nargin > 2
    poi=find(id==codelist);
    if ~isempty(poi)
        id=poi(end);
    else
        id=length(tags)+1;
    end
end

lookup=find(tags==delimiter);
if id<length(lookup)
    tag=tags(lookup(id)+1:lookup(id+1)-1);
end

function id = tag2id(tag,tags,codelist,delimiter)
% function id=tag2id(tag,tags,codelist,delimiter)
%
% Returns the identification code corresponding to a string tag by
% comparison with a list of possible tags given in a string with colon (:)
% separation, another separator can be selected by optional parameter
% delimiter
% if an array codelist is given the identifier is of the same type as the
% elements in codelist, the sequence in the codelist must correspond to the
% one in the tags string
% if codelist is missing, the id is an integer corresponding to the
% position of the tag in string tags
% if tag is not found in tags or if the position is larger than the number
% of elements in codelist, an empty id is given back
%
% tag       (string) tag to be found, example: 'Be'
% tags      (string) tag list, example ':H:He:Li:Be:C:N:O:F:Ne:'
% codelist  [optional] (array of ids), example [1,4,7,9,12,14,16,19,20],
%           can be empty to allow different delimiter without codelist
% delimiter [optional] (char), alternative delimiter, defaults to ':'
%
% id        id selected from array codelist or integer telling the position
%           of tag in tags, for example above:
%           if codelist is missing, id=4
%           if codelist is present, id=9
%
% G. Jeschke, 2009

id=[]; % empty default output
if nargin < 4
    delimiter=':'; % colon as default delimiter
end

etag=[delimiter tag delimiter];
position=strfind(tags,etag);
if position
    id=1+sum(find(tags==delimiter)<position(1));
    if nargin>=3 && ~isempty(codelist)
        if id<=length(codelist)
            id=codelist(id);
        else
            id=[];
        end
    end
end

function show_snake(ofid,pop,sadr,address)
% 
fprintf(ofid,'show %s%s coil %6.3f\n',sadr,address,sqrt(pop));
[chains,residues] = split_address(address);
seg_length = residues(end) - residues(1) + 1;
for c = 1:length(chains)
    first = sprintf('%s(%s)%i',sadr,chains(c),residues(1));
    fprintf(ofid,'colorscheme %s%s sequence %s %i\n',sadr,address,first,seg_length);
end
% for r = 1:length(residues)
%     for c = 1:length(chains)        
%         fprintf(ofid,'color %s(%s)%i %6.3f %6.3f %6.3f\n',sadr,chains(c),residues(r),color_grade(r,length(residues)));
%     end
% end
