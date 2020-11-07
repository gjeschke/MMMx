function [all_files,pop,exceptions] = rd_ensemble_definition(fname)
%
% RD_ENSEMBLE_DEFINITION Read conformer file names and populations
%
%   [all_files,popdef,exceptions] = rd_ensemble_definition(fname)
%   Reads the file names and populations of conformers from an MMMx
%   ensemble definition file
%
% INPUT
% fname     file name, extension .ens is appended if there is none
%
% OUTPUT
% all_files     array of structures with single field (to be consistent
%               with dir)
%               .name file name
% pop           populations
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

[~,~,ext] = fileparts(fname);
if isempty(ext)
    fname = strcat(fname,'.ens');
end

% open the PDB file, return exception if impossible
try
    fid = fopen(fname);
catch exception
    exceptions{1} = exception;
    return
end

if fid == -1
    all_files = [];
    pop = [];
    exceptions = {MException('rd_ensemble_definition:file_could_not_be_opened',...
                                    'File %s could not be opened',fname)};
    return
end

% pre-allocate memory
all_files(10000).name = '';
pop = zeros(1,10000);
C = 0; % initialize number of conformers

while 1
    tline = fgetl(fid);
    if ~ischar(tline) 
        break 
    end
    poi = strfind(tline,'%');
    if ~isempty(poi)
        if poi(1) == 1
            continue
        end
        tline = tline(1:poi-1);
    end
    args = split(tline);
    if length(args) >= 2 % well-formed name/population line
        C = C + 1;
        all_files(C).name = args{1};
        pop(C) = str2double(args{2});
    end
end

fclose(fid);

all_files = all_files(1:C);
pop = pop(1:C);