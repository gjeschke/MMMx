function [pairs,labels,exceptions] = rd_pair_list(fname)
%
% RD_PAIR_LIST Read list of label site pairs (as written by pair_list.m or
%              hetero_pair_list.m
%
%   [pairs,exceptions] = rd_pair_list(fname)
%   Reads site addresses, mean distance, standrad deviation from an MMMx 
%   label site pair file
%
% INPUT
% fname     file name, extension .lst is appended if there is none
%
% OUTPUT
% pairs         array of structures with fields
%               .address1   address of site 1
%               .address2   address of site 2
%               .r_mean     mean distance
%               .r_std      standard deviation of distance
% labels        cell(1,2) string with label (rotamer library) names
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

labels = {};
exceptions{1} = [];
pairs(1000).address1 = '';
pairs(1000).address2 = '';
pairs(1000).r_mean = 0;
pairs(1000).r_std = 0;

% append extension if required
[~,~,ext] = fileparts(fname);
if isempty(ext)
    fname = strcat(fname,'.lst');
end

% open the file, return exception if impossible
try
    fid = fopen(fname);
catch exception
    exceptions{1} = exception;
    pairs = [];
    return
end

if fid == -1
    pairs = [];
    exceptions = {MException('rd_pair_list:file_could_not_be_opened',...
                                    'File %s could not be opened',fname)};
    return
end

P = 0; % initialize number of pairs

while 1
    tline = fgetl(fid);
    if ~ischar(tline) 
        break 
    end
    poi = strfind(tline,'%');
    if ~isempty(poi)
        if poi(1) == 1
            tokpoi = strfind(tline,'label:');
            if ~isempty(tokpoi)
                labels{1} = strtrim(tline(tokpoi+6:end));
                labels{2} = labels{1};
            else
                tokpoi = strfind(tline,'labels:');
                if ~isempty(tokpoi)
                    labels = split(tline(tokpoi+7:end),'|');
                    labels{1} = strtrim(labels{1});
                    labels{2} = strtrim(labels{2});
                end
            end
            continue
        end
        tline = tline(1:poi-1);
    end
    args = split(tline);
    if length(args) >= 4 % well-formed site line
        P = P + 1;
        pairs(P).address1 = args{1};
        pairs(P).address2 = args{2};
        pairs(P).r_mean = str2double(args{3});
        pairs(P).r_std = str2double(args{4});
    end
end

fclose(fid);

pairs = pairs(1:P);