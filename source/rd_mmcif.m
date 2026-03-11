function info = rd_mmcif(fname,missing_comments)
%
% RD_MMCIF Reads mmcif file into Matlab struct
%
%   info = rd_mmcif(fname)
%
% INPUT
% fname     file name, extension .cif is appended if there is none
% missing_comments indicates a non-standrad format without comment lines
%                  for sections, defauts to false
%
% OUTPUT
% info      struct, fields and their length are defined by file content 
%           the current version reads only the 3D structure definition    
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2024: Gunnar Jeschke

info = [];

if ~exist('missing_comments','var')
    missing_comments = false;
end

[~,~,ext] = fileparts(fname);
if isempty(ext)
    fname = strcat(fname,'.cif');
end

% open the file, return exception if impossible
try
    fid = fopen(fname);
catch
    return
end

if fid == -1
    return
end

clear info

tline = fgetl(fid);
if ~missing_comments
    while 1
        if ~ischar(tline)
            break
        end
        poi = strfind(tline,'#');
        if ~isempty(poi)
            sfpoi = 0;
            if poi(1) == 1 % this is a new field for the structure
                tline = fgetl(fid);
                if ~ischar(tline)
                    break
                end
                if length(tline) >= 5 && strcmpi(tline(1:5),'loop_') % this an array field
                    array = true;
                    subfields{100} = '';
                    sfpoi = 0;
                    loop_poi = 0;
                    tline = fgetl(fid);
                    if ~ischar(tline)
                        break
                    end
                else
                    array = false;
                    name_arg = split(tline);
                    if strcmpi(name_arg{1},'_entry.id')
                        info.entry_id = name_arg{2};
                    end
                end
            end
            while tline(1) == '_'
                tline = tline(2:end);
                if array % the only array variable that is currently read are atom sites
                    names = split(tline,'.');
                    cfield = names{1};
                    if strcmpi(cfield,'atom_site')
                        cfield = 'atomsite';
                        subfield = strtrim(names{2});
                        poi = 0;
                        for c = 1:length(subfield)
                            if subfield(c) ~= '_'
                                poi = poi + 1;
                                subfield(poi) = subfield(c);
                            end
                        end
                        subfield = subfield(1:poi);
                        info.(cfield)(1).(subfield) = '';
                        sfpoi = sfpoi + 1;
                        subfields{sfpoi} = subfield;
                    else
                        tline = fgetl(fid);
                        if ~ischar(tline)
                            break
                        end
                        continue
                    end
                else
                    tline = fgetl(fid);
                    if ~ischar(tline)
                        break
                    end
                    continue
                    % cfield = name_arg{1};
                    % arg = name_arg{2};
                    % info.(cfield) = arg;
                end
                tline = fgetl(fid);
                if ~ischar(tline)
                    break
                end
            end
            if sfpoi == 0
                continue
            end
            while tline(1) ~= '#'
                loop_poi = loop_poi + 1;
                args = split(tline);
                for sf = 1:sfpoi
                    subfield = subfields{sf};
                    info.(cfield)(loop_poi).(subfield) = args{sf};
                end
                tline = fgetl(fid);
                if ~ischar(tline)
                    break
                end
            end
        else
            tline = fgetl(fid);
            if ~ischar(tline)
                break
            end
        end
    end
else
    while 1
        if ~ischar(tline)
            break
        end
        if length(tline) >= 5 && strcmpi(tline(1:5),'loop_') % this an array field
            array = true;
            subfields{100} = '';
            sfpoi = 0;
            loop_poi = 0;
            tline = fgetl(fid);
            if ~ischar(tline)
                break
            end
        else
            array = false;
            name_arg = split(tline);
            if strcmpi(name_arg{1},'_entry.id')
                info.entry_id = name_arg{2};
            end
        end
        if array
            while tline(1) == '_'
                tline = tline(2:end);
                if array % the only array variable that is currently read are atom sites
                    names = split(tline,'.');
                    cfield = names{1};
                    if strcmpi(cfield,'atom_site')
                        cfield = 'atomsite';
                        subfield = strtrim(names{2});
                        poi = 0;
                        for c = 1:length(subfield)
                            if subfield(c) ~= '_'
                                poi = poi + 1;
                                subfield(poi) = subfield(c);
                            end
                        end
                        subfield = subfield(1:poi);
                        info.(cfield)(1).(subfield) = '';
                        sfpoi = sfpoi + 1;
                        subfields{sfpoi} = subfield;
                    else
                        tline = fgetl(fid);
                        if ~ischar(tline)
                            break
                        end
                        continue
                    end
                else
                    tline = fgetl(fid);
                    if ~ischar(tline)
                        break
                    end
                    continue
                    % cfield = name_arg{1};
                    % arg = name_arg{2};
                    % info.(cfield) = arg;
                end
                tline = fgetl(fid);
                if ~ischar(tline)
                    break
                end
            end
            if sfpoi == 0
                continue
            end
            while length(tline) > 5 && tline(1) ~= '#'
                loop_poi = loop_poi + 1;
                args = split(tline);
                for sf = 1:sfpoi
                    subfield = subfields{sf};
                    info.(cfield)(loop_poi).(subfield) = args{sf};
                end
                tline = fgetl(fid);
                if ~ischar(tline)
                    break
                end
            end
        else
            tline = fgetl(fid);
        end
    end
end

fclose(fid);

