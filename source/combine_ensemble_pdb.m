function exceptions = combine_ensemble_pdb(description,fname)
%
% COMBINE_ENSEMBLE_PDB Combine several PDB files to an ensemble file
%
%   exceptions = COMBINE_ENSEMBLE_PDB(decription,fname)
%   Writes an ensemble PDB file with population information, the ensemble
%   is described by an MMMx .ens file, as written by module_ensemble_fit
%   the individual PDB files must be on the Matlab path
%
% INPUT
% description  file name of the ensemble description
% fname        file name for output PDB file, default generated from
%              descriptor file
%
% OUTPUT
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
exceptions = {[]};

[all_files,pop,exceptions_rd_def] = rd_ensemble_definition(description);

if ~isempty(exceptions_rd_def{1})
    exceptions = exceptions_rd_def;
    return
end

if isempty(all_files) % error, if ensemble is empty
    exceptions = {MException('combine_ensemble_pdb:empty_ensemble',...
                                    'Ensemble list %s is empty',description)};
                                return
end

if length(all_files) == 1 % error, if ensemble has only one conformer
    exceptions = {MException('combine_ensemble_pdb:no_ensemble',...
                                    'Ensemble list %s has only one conformer',description)};
                                return
end

% make output file name, if none

if ~exist('fname','var') || isempty(fname)
    poi = strfind(description,'.ens');
    if ~isempty(poi)
        fname = [description(1:poi-1) '_ensemble.pdb'];
    else
        fname = [decsription '_ensemble.pdb'];
    end
end

if ~contains(fname,'.pdb')
    fname = [fname '.pdb'];
end

% try to open output file
% open the PDB file, return exception if impossible
try
    ofid = fopen(fname,'wt');
catch exception
    exceptions{1} = exception;
    return
end

% open the PDB file, return exception if impossible
try
    fid = fopen(all_files(1).name);
catch exception
    exceptions{1} = exception;
    return
end

if fid == -1
    exceptions = {MException('combine_ensemble_pdb:conformer_missing',...
                                    'Conformer file %s could not be opened',all_files(1).name)};
    return
end

found_header = false;
atom_section = false;
model_opened = false;
pop_written = false;
while 1
    tline = fgetl(fid);
    if ~ischar(tline) 
        break 
    end
    if strcmp(tline(1:6),'HEADER') % check whether file has a header line
        found_header = true;
    end
    % check wheter ATOM section has started
    if strcmp(tline(1:4),'ATOM') % check whether we have encountered an atom
        atom_section = true;
    end    
    % check whether REMARK section has ended
    if strcmp(tline(1:6),'SEQRES') || strcmp(tline(1:5),'SHEET') || strcmp(tline(1:5),'HELIX')...
       || atom_section
        % write HEADER if there is none
        if ~found_header
            fprintf(ofid,'HEADER    MMMX MINIMAL ENSEMBLE PDB FILE          %s   %s\n', ...
                datestr(datetime,'dd-mmm-yy'),pdbid);
        end
        if ~pop_written
            % write conformer population information if requested
            pdbline = 'REMARK 400  POPULATIONS';
            fprintf(ofid,'%s\n',pad(pdbline,80));
            for c = 1:length(pop)
                pdbline = sprintf('REMARK 400   MODEL %9i POPULATION %8.4f',c,pop(c));
                fprintf(ofid,'%s\n',pad(pdbline,80));
            end
            pop_written = true;
        end
    end
    if ~atom_section
        if ~strcmpi(tline(1:6),'END   ') && ~strcmpi(tline(1:6),'MASTER')
            fprintf(ofid,'%s\n',tline);
        end
    end
    if atom_section 
        if ~model_opened
            pdbline = sprintf('MODEL %8i',1);
            fprintf(ofid,'%s\n',pad(pdbline,80));
            model_opened = true;
        end
        if ~strcmpi(tline(1:6),'END   ') && ~strcmpi(tline(1:6),'MASTER')
            fprintf(ofid,'%s\n',tline);
        end
    end
end
fclose(fid);

pdbline = 'ENDMDL';
fprintf(ofid,'%s\n',pad(pdbline,80));

for c = 2:length(all_files) % all conformers
    % open the PDB file, return exception if impossible
    try
        fid = fopen(all_files(c).name);
    catch exception
        exceptions{1} = exception;
        return
    end
    
    if fid == -1
        exceptions = {MException('combine_ensemble_pdb:conformer_missing',...
            'Conformer file %s could not be opened',all_files(c).name)};
        return
    end
    
    atom_section = false;
    model_opened = false;
    while 1
        tline = fgetl(fid);
        if ~ischar(tline)
            break
        end
        % check wheter ATOM section has started
        if strcmp(tline(1:4),'ATOM') % check wehther file has a header line
            atom_section = true;
        end
        if atom_section
            if ~model_opened
                pdbline = sprintf('MODEL %8i',c);
                fprintf(ofid,'%s\n',pad(pdbline,80));
                model_opened = true;
            end
            if ~strcmpi(tline(1:6),'END   ') && ~strcmpi(tline(1:6),'MASTER')
                fprintf(ofid,'%s\n',tline);
            end
        end
    end
    fclose(fid);
    pdbline = 'ENDMDL';
    fprintf(ofid,'%s\n',pad(pdbline,80));
end

pdbline = 'END';
fprintf(ofid,'%s\n',pad(pdbline,80));
fclose(ofid);