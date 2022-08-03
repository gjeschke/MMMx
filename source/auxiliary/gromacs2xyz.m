function gromacs2xyz(fname,options)
% GROMACS2XYZ Converts a GROMACS coordinate file to an xyz or Chem3D cc1
% file, input is in nm, output is in Angstroem, output velocities are in
% Angstroem/ps if requested
%
% gromacs2xyz(options) 
%
% fname     input filename, extension .gro is appended, if there is none 
% options   structure of conversion options with fields
%           .cc1    Chem3D .cc1 format, where the title line is skipped
%                   Boolean flag, optional, defaults to false
%           .fname  filename for output, string, optional, defaults to
%                   input filename with extension .xyz or .cc1
%           .sort   cell string of molecule types (three-letter codes or 
%                   five-letter names),
%                   output is sorted in the order given and only molecule
%                   types in this list are output, optional, default is
%                   output of all atoms in the same order as in the input
%                   file
%           .title  string, optional, title line for output, defaults to
%                   title in input file, is neglected if output format is
%                   .cc1
%           .v      Boolean flag that requests velocity output, defaults to
%                   false, if the input file does not specify velocities,
%                   they are set to NaN in the output
%           .assign element assignment for atom types, cell(n,2) of
%                   strings, where the first string of each row is an atom
%                   type and the second string an element symbol, if
%                   missing, only elements with single-letter symbols are
%                   correctly recognized, if present, all atom types in the
%                   input file must be included
%
% GROMACS format specification taken from https://manual.gromacs.org/archive/5.0.3/online/gro.html
%
% GROMACS atom names do not allow for unambiguous assignment of the element
% without access to the atomtypes.atp file of the force field (and even
% with this file, element assignment could fail for atom groups), if in
% doubt, specify an assignment list in options

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2022: Gunnar Jeschke

% analyze input and set defaults
[filepath,bname,ext] = fileparts(fname);
if isempty(ext)
    ext = '.gro';
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

if isfield(options,'assign') && ~isempty(options.assign)
    assign = true;
    [types,~] = size(options.assign);
else
    assign = false;
end

% coordinate input section, coordinates are stored either sequentially or
% per molecule type, depending on whether output should be sorted

fid = fopen(fname,'rt');
if fid == -1
    error('gromacs2xyz:no_file','Input file %s could not be opened\n',fname);
end

if isempty(title)
    title = fgetl(fid);
else
    fgetl(fid);
end
    
n_atoms = str2double(fgetl(fid));

if ~sorted
    xyz = zeros(n_atoms,3);
    v = zeros(n_atoms,3);
    elements = cell(n_atoms,1);
    atoms = n_atoms;
else
    coordinates.atoms = n_atoms;
end

fprintf(1,'Reading %i atom coordinates\n',n_atoms);
for k = 1:n_atoms
    % read and analyze one input line
    tline = fgetl(fid);
    resnum = str2double(tline(1:5));
    restype = strtrim(tline(6:10)); 
    atomtype = strtrim(tline(11:15));
    atnum = str2double(tline(16:20));
    x = 10*str2double(tline(21:28));
    y = 10*str2double(tline(29:36));
    z = 10*str2double(tline(37:44));
    if length(tline) >= 68
        vx = 10*str2double(tline(45:52));
        vy = 10*str2double(tline(53:60));
        vz = 10*str2double(tline(61:68));
    else
        vx = NaN;
        vy = NaN;
        vz = NaN;
    end
    % determine element
    if assign
        element = '';
        t = 0;
        while isempty(element) && t < types
            t = t + 1;
            if strcmpi(options.assign(t,1),atomtype)
                element = options.assign(t,2);
            end
        end
        if isempty(element)
            fclose(fid);
            error('gromacs2xyz:unknown_atom_type','Atom type %s not specified in assignment list\n',atomtype);
        end
    else
        element = atomtype(1);
    end
    if sorted
        % initialize new residue type if required
        if ~isfield(coordinates,restype)
            coordinates.(restype).N = 0;
            coordinates.(restype).xyz = zeros(n_atoms,3);
            coordinates.(restype).v = zeros(n_atoms,3);
            coordinates.(restype).elements = cell(n_atoms,1);
            coordinates.(restype).atnum = zeros(n_atoms,1);
            coordinates.(restype).resnum = zeros(n_atoms,1);
        end
        coordinates.(restype).N = coordinates.(restype).N + 1;
        n = coordinates.(restype).N;
        coordinates.(restype).xyz(n,:) = [x,y,z]; 
        coordinates.(restype).v(n,:) = [vx,vy,vz]; 
        coordinates.(restype).elements{n} = element; 
        coordinates.(restype).atnum(n) = atnum; 
        coordinates.(restype).resnum(n) = resnum; 
    else
        xyz(k,:) = [x,y,z]; 
        v(k,:) = [vx,vy,vz]; 
        elements{k} = element;
    end
end

fclose(fid);

% determine number of output atoms if output is sorted

if sorted
    atoms = 0;
    for t = 1:length(options.sort)
        if isfield(coordinates,options.sort{t})
            fprintf(1,'%i atom coordinates of residue type %s will be written\n',coordinates.(options.sort{t}).N,options.sort{t});
            atoms = atoms + coordinates.(options.sort{t}).N;
        else
            error('gromacs2xyz:missing_residue','Residue type %s required in output, but not present in input\n',options.sort{t});
        end
    end
end

% write output

fid = fopen(options.fname,'wt');
if fid == -1
    error('gromacs2xyz:no_output','Output file %s could not be opened\n',options.fname);
end

if ~options.cc1
    fprintf(fid,'%s\n',title);
end

fprintf(fid,'%i\n',atoms);

if sorted
    for t = 1:length(options.sort)
        restype = options.sort{t};
        for a = 1:coordinates.(restype).N
            fprintf(fid,'%3s',coordinates.(restype).elements{a});
            fprintf(fid,'%12.6f',coordinates.(restype).xyz(a,:));
            if options.v
                fprintf(fid,'%12.6f',coordinates.(restype).v(a,:));
            end
            fprintf(fid,'\n');
        end
    end
else
    for a = 1:atoms
        fprintf(fid,'%3s',elements{a});
        fprintf(fid,'%12.6f',xyz(a,:));
        if options.v
            fprintf(fid,'%12.6f',v(a,:));
        end
        fprintf(fid,'\n');
    end
end

fclose(fid);