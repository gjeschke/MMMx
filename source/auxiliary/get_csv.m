function [data,description] = get_csv(fname,format)
%
% GET_CSV    Read data from a comma-separated value file
%
%   data = GET_CSV(fname,format)
%   Reads data from a comma-separated value file. The data may or may not
%   contain a title line and may or may not contain a unit line. The body
%   lines can be mixed format or purely numerical. The number of items per
%   row must be the same in any row as in the first row.
%
% INPUT
% fname         name of file, extension .csv is appended if there is none
% format        string, specifies input and output format in the form
%               <input>.<output>
%               <input>     data        only data, no header
%                           header      header line with variable names only
%                           units       header line with variable names and
%                                       units line
%               <output>    table       Matlab table as output
%                           cell        cell array as output
%                           numerical   numerical array
%               Example:    header.numerical
%               Default:    header.table            
%
% OUTPUT
% data          data in the requested output format, table format includes
%               the description
% description   either cell (1,m) of column names or cell(2,m) with column
%               names in the first row and units in the second
%               row
%
% Remarks:
%
% This function reads data output by put_csv.m as well as typical .csv data
%
% In output format 'numerical', output is empty if any body variable is not
% a number
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

data = [];
description = {''};

[pname,bname,ext] = fileparts(fname); 
if isempty(ext)
    fname = fullfile(pname,strcat(bname,'.csv'));
end

if ~exist('format','var')
    format = 'header.table';
end

args = split(format,'.');
input_type = args{1};
output_type = args{2};

% read data as a Matlab table
fid = fopen(fname,'rt');

nl = 0;
rows = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    nl = nl + 1;
    rows = rows + 1;
    args = split(tline,',');
    fprintf(1,'Line %i has %i arguments\n',nl,length(args));
    types = cell(1,length(args));
    for k = 1:length(args)
        types{k} = 'cellstr';
    end
    if nl == 1
        data = table('Size',[10000,length(args)],'VariableTypes',types);
        if strcmpi(input_type,'header') || strcmpi(input_type,'units')
            data.Properties.VariableNames = args'; 
            description = args';
            rows = rows - 1;
        end
    end
    if nl == 2
        if strcmpi(input_type,'units')
            data.Properties.VariableUnits = args';
            description = [description;args'];
            rows = rows - 1;
        end
    end
    if rows > 0
        data(rows,:) = args'; %#ok<AGROW>
    end
end

        
fclose(fid);

data = data(1:rows,:);
if strcmpi(output_type,'cell')
    data = data.Variables;
    return
end

% convert numerical data cells to double 

[m,n] = size(data);
numerical = ones(1,n);
for row = 1:m
    for column = 1:n
        trial = str2double(data{row,column});
        if isnan(trial)
            numerical(column) = 0;
        else
            data{row,column} = {trial};
        end
    end
end

if sum(numerical) > 0
    types = cell(1,n);
    for column = 1:n
        if numerical(column)
            types{column} = 'double';
        else
            types{column} = 'cellstr';
        end
    end
    newdata = table('Size',[m,n],'VariableTypes',types);
    newdata.Properties.VariableNames = data.Properties.VariableNames;
    newdata.Properties.VariableUnits = data.Properties.VariableUnits;
    for row = 1:m
        for column = 1:n
            if numerical(column)
                item = data{row,column};
                newdata{row,column} = item{1};
            else
                newdata{row,column} = data{row,column};
            end
        end
    end
    data = newdata;
end

if strcmpi(output_type,'numerical')
    if sum(numerical) == length(numerical)
        data = data.Variables;
    else
        data = [];
    end
end

