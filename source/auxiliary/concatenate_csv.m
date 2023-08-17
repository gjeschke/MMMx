function concatenate_csv(fname,filelist)
%
% CONCATENATE_CSV    Concatenates several comma-separated value files
%
%   CONCATENATE_CSV(fname,filelist)
%   Concatenates several comma-separated value files into a single file in
%   order to conform to the (silly) source-data format for figures in
%   Nature group journals
%
% INPUT
% fname         name of file, extension .csv is appended if there is none
% filelist      list of file names for files to be included, cellstring
%
% Remarks:
%
% each section of an input file is preceded by a CSV comment line (starting
% with #) that contains only the original file name
%
% if the name in the file list does not have an extension, .csv is assumed
%
% any text file can be an input file, input files are copied verbatim into
% the output file
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

[pname,bname,ext] = fileparts(fname); 
if isempty(ext)
    fname = fullfile(pname,strcat(bname,'.csv'));
end

ofid = fopen(fname,'wt');

for k = 1:length(filelist)
 
    iname = filelist{k};
    [pname,bname,ext] = fileparts(iname);
    if isempty(ext)
        iname = fullfile(pname,strcat(bname,'.csv'));
    end

    fprintf(ofid,'# %s\n',iname);
    % read/write data for a single file
    ifid = fopen(iname,'rt');
    
    while 1
        tline = fgetl(ifid);
        if ~ischar(tline), break, end
        fprintf(ofid,'%s\n',tline);
    end
    
    
    fclose(ifid);

end

fclose(ofid);