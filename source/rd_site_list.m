function [sites,label,chains,exceptions] = rd_site_list(fname)
%
% RD_SITE_LIST Read list of label sites (as written by site_scan.m)
%
%   [sites,exceptions] = rd_site_list(fname)
%   Reads site addresses, number of rotamers, partition function, label
%   position rmsd and original residue from an MMMx label site file
%
% INPUT
% fname     file name, extension .lst is appended if there is none
%
% OUTPUT
% sites         array of structures with fields
%               .address    site address
%               .rotamers   number of rotamers
%               .Z          partition function
%               .rmsd       label position rmsd (Angstroem)
%               .residue    name of original residue
% label         name of label (rotamer library)
% chains        string with single-letter chain identifiers, can be '*' for
%               all chains
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

label = '';
exceptions{1} = [];
sites(1000).address = '';
sites(1000).rotamers = 0;
sites(1000).Z = 0;
sites(1000).rmsd = 0;
sites(1000).residue = '';

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
    sites = [];
    return
end

if fid == -1
    sites = [];
    exceptions = {MException('rd_site_list:file_could_not_be_opened',...
                                    'File %s could not be opened',fname)};
    return
end

S = 0; % initialize number of sites

tline = fgetl(fid);
pa = strfind(tline,'(');
pe = strfind(tline,')');
if ~isempty(pa) && ~isempty(pe)
    chains = tline(pa+1:pe-1);
else
    chains = '*';
end
poi = strfind(tline,'%');
if ~isempty(poi)
    if poi(1) == 1
        tokpoi = strfind(tline,'label:');
        if ~isempty(tokpoi)
            label = strtrim(tline(tokpoi+6:end));
        end
    end
end

while 1
    tline = fgetl(fid);
    if ~ischar(tline) 
        break 
    end
    poi = strfind(tline,'%');
    if ~isempty(poi)
        tline = tline(1:poi-1);
    end
    args = split(tline);
    if length(args) >= 5 % well-formed site line
        S = S + 1;
        sites(S).address = args{1};
        sites(S).rotamers = str2double(args{2});
        sites(S).Z = str2double(args{3});
        sites(S).rmsd = str2double(args{4});
        name = args{5};
        sites(S).residue = name(2:end-1);
    end
end

fclose(fid);

sites = sites(1:S);