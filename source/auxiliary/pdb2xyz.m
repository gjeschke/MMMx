function pdb2xyz(fname,options)
% PDB2XYZ Converts a PDB to an xyz or Chem3D cc1
% file
%
% pdb2xyz(options) 
%
% This version reads coordinates of only H, O, C, and N atoms
%
% fname     input filename, extension .pdb is appended, if there is none 
% options   structure of conversion options with fields
%           .cc1    Chem3D .cc1 format, where the title line is skipped
%                   Boolean flag, optional, defaults to false
%           .fname  filename for output, string, optional, defaults to
%                   input filename with extension .xyz or .cc1
%           .title  string, optional, title line for output, defaults to
%                   title in input file, is neglected if output format is
%                   .cc1
%           .v      Boolean flag that requests velocity output, defaults to
%                   false, if the input file does not specify velocities,
%                   they are set to NaN in the output
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2022: Gunnar Jeschke

% analyze input and set defaults
[filepath,bname,ext] = fileparts(fname);
if isempty(ext)
    ext = '.pdb';
end
fname = fullfile(filepath,strcat(bname,ext));

if ~exist('options','var') || ~isfield(options,'cc1') || isempty(options.cc1)
    options.cc1 = false;
end

if ~isfield(options,'fname') || isempty(options.fname)
    if options.cc1
        options.fname = fullfile(filepath,strcat(bname,'.cc1'));
    else
        options.fname = fullfile(filepath,strcat(bname,'.xyz'));
    end
end

if isfield(options,'sort') && ~isempty(options.sort)
    sorted = true;
    for s = 1:length(options.sort)
        options.sort{s} = strtrim(options.sort{s});
    end
else
    sorted = false;
end

if isfield(options,'title') && ~isempty(options.title)
    title = options.title;
else
    title = '';
end

if ~isfield(options,'v') || isempty(options.v)
    options.v = false;
end

% coordinate input section, coordinates are stored either sequentially or
% per molecule type, depending on whether output should be sorted

fid = fopen(fname,'rt');
if fid == -1
    error('pdb2xyz:no_file','Input file %s could not be opened\n',fname);
end

if isempty(title)
    title = 'From PDB file read by MMMx';
end

n_atoms = 50000;
xyz = zeros(n_atoms,3);
elements = cell(n_atoms,1);

fprintf(1,'Reading atom coordinates\n');
n_atoms = 0;
while 1
    % read and analyze one input line
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if length(tline) < 4
        continue;
    end
    if ~strcmpi(tline(1:4),'ATOM') && length(tline) < 6
        continue
    elseif ~strcmpi(tline(1:4),'ATOM') && ~strcmpi(tline(1:6),'HETATM')
        continue
    end
    x = str2double(tline(31:38));
    y = str2double(tline(39:46));
    z = str2double(tline(47:54));
    % determine element
    element = tline(78);
    if element == 'H' || element == 'O' || element == 'N' || element == 'C'
        n_atoms = n_atoms + 1;
    else
        continue;
    end
    xyz(n_atoms,:) = [x,y,z];
    elements{n_atoms} = element;
end

fclose(fid);

% write output

fid = fopen(options.fname,'wt');
if fid == -1
    error('pdb2xyz:no_output','Output file %s could not be opened\n',options.fname);
end

if ~options.cc1
    fprintf(fid,'%s\n',title);
end

fprintf(fid,'%i\n',n_atoms);

for a = 1:n_atoms
    fprintf(fid,'%3s',elements{a});
    fprintf(fid,'%12.6f',xyz(a,:));
    fprintf(fid,'\n');
end


fclose(fid);