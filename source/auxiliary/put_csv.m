function put_csv(fname,data,description,significant)
%
% PUT_CSV    Write data to a comma-separated value file
%
%   PUT_CSV(fname,data,description,significant)
%   Write a data set with rows of a uniform format to a comma-separated
%   value file. Data can be a numerical array or a cell array. Format can
%   be specified or is determined from data. The first line of the file
%   specifies contents of columns and is derived from input description.
%
% INPUT
% fname         name of file, extension .csv is appended if there is none
% data          either numerical matrix (n,m) or cell array (n,m) or table
% description   either cell (1,m) of column names or cell(2,m) with column
%               names in the first row and format specifiers in the second
%               row, for data in table format, description is overwritten
% significant   number of significant digits for floats, defaults to 6
%
% Remarks:
%
% Automatic formatting recognizes integers and whether a real number can be
% written with a reasonable number of characters in decimal-point format or
% needs to be written in scientific notation 
%
% For data in Matlab table format, a second headline row is written if
% variable units are specified

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

if ~exist('significant','var') || isempty(significant)
    significant = 6;
end

[pname,bname,ext] = fileparts(fname); 
if isempty(ext)
    fname = fullfile(pname,strcat(bname,'.csv'));
end

if istable(data)
    description = data.Properties.VariableNames;
end

header = {};
format = {};
if exist('description','var') && ~isempty(description)
    [n,m] = size(description);
    header = description(1,:);
    if n > 1
        format = description(2,:);
    end
end

fid = fopen(fname,'wt');
% write headline, if any
if ~isempty(header)
    fprintf(fid,'%s',header{1});
    for column = 2:m
        fprintf(fid,',%s',header{column});
    end
    fprintf(fid,'\n');
end

% write units, if data is in table format and units are specified
if istable(data) && ~isempty(data.Properties.VariableUnits)
    fprintf(fid,'%s',data.Properties.VariableUnits{1});
    for column = 2:m
        fprintf(fid,',%s',data.Properties.VariableUnits{column});
    end
    fprintf(fid,'\n');    
end

% write data
[n,m1] = size(data);
if exist('m','var') && m1 ~= m
    fprintf(2,'Warning: Headline has %i columns, but data has %i columns\n',m,m1);
end
for row = 1:n
    if ~isempty(format)
        if istable(data)
            item = data{row,1};
            if iscell(item)
                item = item{1};
            end
        elseif iscell(data)
            item = data{row,1};
        else
            item = data(row,1);
        end
        fprintf(fid,format{1},item);
        for column = 2:m1
            if istable(data)
                item = data{row,column};
                if iscell(item)
                    item = item{1};
                end
            elseif iscell(data)
                item = data{row,column};
            else
                item = data(row,column);
            end
            formstr = sprintf(',%s',format{column});
            fprintf(fid,formstr,item);
        end
    else
        if istable(data)
            item = data{row,1};
            if iscell(item)
                item = item{1};
            end
        elseif iscell(data)
            item = data{row,1};
        else
            item = data(row,1);
        end
        fprintf(fid,'%s',autoformat(item,significant));
        for column = 2:m1
            if istable(data)
                item = data{row,column};
                if iscell(item)
                    item = item{1};
                end
            elseif iscell(data)
                item = data{row,column};
            else
                item = data(row,column);
            end
            fprintf(fid,',%s',autoformat(item,significant));
        end
    end
    fprintf(fid,'\n');
end
fclose(fid);

function prt_string = autoformat(item,significant)

prt_string = '###ERROR###';
if iscell(item)
    item = item{1};
end
[m,n] = size(item);
if ~ischar(item)
    if m ~=1 || n~= 1
        return
    end
end
if isinteger(item)
    prt_string = sprintf('%i',item);
elseif isfloat(item)
    format = sprintf('%%%i.%ig',8+significant,significant);
    prt_string = sprintf(format,item);
elseif ischar(item)
    prt_string = item;
end
prt_string = strtrim(prt_string);