function extract_from_gromacs_trr(fname,template,frames,options)
% EXTRACT_FROM_GROMACS_TRR   Extracts coordinates from the text dump of a
% GROMACS .trr binary trajectory file
%
% use something like 
% gmx dump -f wg_box_50.trr > wg_box_50_trj.dat
% to create the text dump
%
% extract_from_gromacs_trr(fname,template,frames,options) 
%
% fname     input filename, extension .gro is appended, if there is none
% template  GROMACS .gro file with information on the system
% frames    vector of numbers of frames that are to be extracted
% options   structure of conversion options with fields
%           .cc1    Chem3D .cc1 format, where the title line is skipped
%                   Boolean flag, optional, defaults to false
%           .fname  filename for output, string, optional, defaults to
%                   input filename with extension .xyz or .cc1 , the frame
%                   numbers are included automatically, one file is created
%                   for each frame
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

[filepath,bname] = fileparts(fname);


if ~exist('options','var') || ~isfield(options,'cc1') || isempty(options.cc1)
    options.cc1 = false;
end

if ~isfield(options,'fname') || isempty(options.fname)
    options.fname = fullfile(filepath,bname);
end

if isfield(options,'title') && ~isempty(options.title)
    title = options.title;
else
    title = '';
end

[sorting,elements,~,n_atoms] = gromacs2xyz(template,options);


fid = fopen(fname,'rt');
if fid == -1
    error('extract_from_gromacs_trr:no_file','Input file %s could not be opened\n',fname);
end

at_frame = 0;
for f = 1:length(frames)
    at_frame = at_frame + 1;
    while at_frame < frames(f) % skip to next frame
        skip = 8 + 2*n_atoms; % number of lines to be skipped per frame
        while skip > 0
            fgetl(fid);
            skip = skip - 1;
        end
    end
    for l = 1:7 % skip to beginning of coordinate section
        fgetl(fid);
    end
    fprintf(1,'Reading %i atom coordinates of frame %i\n',n_atoms,frames(f));
    fcoor = zeros(n_atoms,3);
    for k = 1:n_atoms
        % read and analyze one input line
        tline = fgetl(fid);
        ccoor = 10*str2num(tline(17:56)); %#ok<ST2NM>
        fcoor(k,:) = ccoor;
    end
    
    % determine number of output atoms if output is sorted
    
    fcoor = fcoor(sorting,:);
    % write output
    
    myname = sprintf('%s_frame_%i',options.fname,frames(f));
    if options.cc1
        myname = strcat(myname,'cc1');
    else
        myname = strcat(myname,'xyz');
    end
    ofid = fopen(myname,'wt');
    if ofid == -1
        error('extract_from_gromacs_trr:no_output','Output file %s could not be opened\n',myname);
    end
    
    if ~options.cc1
        fprintf(ofid,'%s\n',title);
    end
    
    fprintf(ofid,'%i\n',length(sorting));
    for a = 1:length(sorting)
        fprintf(ofid,'%3s',elements{a});
        fprintf(ofid,'%12.6f',fcoor(a,:));
        %     if options.v
        %         fprintf(fid,'%12.6f',v(a,:));
        %     end
        fprintf(ofid,'\n');
    end
    fclose(ofid);
end
fclose(fid);

