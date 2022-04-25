function [entity,exceptions] = get_cyana_pdb(fname,options,entity)
%
% GET_CYANA_PDB Load structure from CYANA dialect PDB file
%
%   [entity,exceptions] = GET_CYANA_PDB(fname)
%   Returns an (entity) structure in MMMx:atomic representation
%
%   [entity,exceptions] = GET_CYANA_PDB(fname,options)
%   Reads topology and Cartesian coordinates of a biological entity
%   from a CYANA dialect PDB file
%
%   [entity,exceptions] = GET_CYANA_PDB(fname,options,entity)
%   Adds Cartesian coordinates of a conformer to an input biological entity
%   by reading from a CYANA dialect PDB file
%
% INPUT
% ident     PDB identifier, if four characters without extension, otherwise
%           file name
% options   structure with requests for extended information, flags
%           .dssp   DSSP (Kabsch/Sanders) secondary structure
%                   defaults to false, DSSP is not performed if there were
%                   insertion codes or non-positive residue numbers
%           .name   optional name for the entity, defaults: PDB identifier
%                   upon download or if header line exists; MMMx otherwise
%           .maxch  maximum number of chains
% entity    optional, if present, models from the PDB file are added as
%           conformers to an existing entity, the caller is responsible for
%           consistency of primary structure of the conformers
%
% OUTPUT
% entity       entity structure in MMMx:atomic representation
% exceptions   error message if something went wrong, entity is empty for
%              errors but not for warnings, for warnings, only the last one
%              is reported, cell containing an empty array, if no exception
%
% Realization:
%
% -     the reader skips residue types that are introduced purely for
%       restraint-based optimization purposes:
%       PL, NL, ML, LL, LL2, LLL, LL5, LP, LN, LM,ION, ORI
% -     the reader converts modified CYS residues to CYS:
%       CYSS, CYSM
% -     the reader converts non-standard nucleotide rseidue names to
%       standard names:
%       RADE to A
%       RCYT to C
%       RGUA to G
%       URA  to U
%       ADE to DA
%       CYT to DC
%       GUA to DG
%       THY to DT


% This file is a part of MMMx. License is MIT (see LICENSE.md). 
%
% G. Jeschke, 2021-2022

exceptions{1} = [];
warnings = 0; % counter for warnings

if ~exist('options','var')
    options = [];
end

if ~exist('entity','var')
    entity = [];
end

maxchain = 26;
if isfield(options,'maxch') && ~isempty(options.maxch)
    maxchain = options.maxch;
end

if ~contains(fname,'.')
    ifid = fopen(strcat(fname,'.pdb'));
else
    ifid = fopen(fname);
end
if ifid == -1
    warnings = warnings + 1;
    exceptions{warnings} = MException('get_cyana_pdb:input_file_missing',...
        'input file %s could not be opened',fname);
    return;
end

tmpname = strcat(fname,'_corr.pdb');
ofid=fopen(tmpname,'wt');
if ofid==-1
    warnings = warnings + 1;
    exceptions{warnings} = MException('get_cyana_pdb:file_not_writeable',...
        'output PDB file %s could not be written',tmpname);
    return;
end

chain_tags = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ';
cres = -1000;
chain = 0;
atnum = 0;
is_atom_line = false;
while 1 && chain <= maxchain
    tline = fgetl(ifid);
    if ~ischar(tline), break, end
    if length(tline)<26 % catches too short lines
        if length(tline) >= 6
            record = tline(1:6);
            if strcmp(record,'MODEL ')
                cres = -1000;
                chain = 0;
                atnum = 0;
            end
        end
        fprintf(ofid,'%s\n',tline);
    else
        record = tline(1:6);
        restype = tline(18:21);
        switch restype
            case ' PL '
                skipped = true;
            case ' LN '
                skipped = true;
            case ' ML '
                skipped = true;
            case ' NL '
                skipped = true;
            case 'LL2 '
                skipped = true;
            case 'LL5 '
                skipped = true;
            case 'LLL '
                skipped = true;
            case ' LP '
                skipped = true;
            case ' LL '
                skipped = true;
            case ' LM '
                skipped = true;
            case 'ION '
                skipped = true;
            case 'ORI '
                skipped = true;
            otherwise
                skipped = false;
        end
        if strcmp(record,'MODEL ')
            cres = -1000;
            chain = 0;
            atnum = 0;
        end
        if strcmp(record,'ATOM  ') || strcmp(record,'HETATM')
            is_atom_line = true;
            resnum = str2double(tline(23:26));
            if resnum ~= cres && resnum ~= cres + 1 && ~skipped
                if chain > 0 % insert chain termination record
                    atnum = atnum + 1;
                    tlinex = tline0;
                    tlinex(1:6) = 'TER   ';
                    tlinex(7:11) = sprintf('%5i',atnum);
                    tlinex = tlinex(1:27);
                    tlinex(13:16) = '    ';
                    tlinex(22) = chain_tags(chain);
                    fprintf(ofid,'%s\n',tlinex);
                end
                chain = chain + 1; % next chain
            end
            if chain > maxchain
                continue
            end
            atnum = atnum + 1;
            tline(7:11) = sprintf('%5i',atnum);
            tline(22) = chain_tags(chain);
            if strcmp(record,'HETATM')
                cres = resnum;
                fprintf(ofid,'%s\n',tline);
            else % wrong ATOM records still need to be fixed
                restype = tline(18:21);
                switch restype
                    case ' PL '
                        atnum = atnum - 1;
                    case ' ML '
                        atnum = atnum -1;
                    case ' NL '
                        atnum = atnum -1;
                    case ' LL '
                        atnum = atnum -1;
                    case 'LL2 '
                        atnum = atnum - 1;
                    case 'LLL '
                        atnum = atnum - 1;
                    case 'LL5 '
                        atnum = atnum - 1;
                    case ' LN '
                        atnum = atnum -1;
                    case ' LP '
                        atnum = atnum -1;
                    case ' LM '
                        atnum = atnum -1;
                    case 'ION '
                        atnum = atnum -1;
                    case 'ORI '
                        atnum = atnum -1;
                    case 'RGUA'
                        tline(18:21) = '  G ';
                        cres = resnum;
                        fprintf(ofid,'%s\n',tline);
                    case 'RADE'
                        tline(18:21) = '  A ';
                        cres = resnum;
                        fprintf(ofid,'%s\n',tline);
                    case 'URA '
                        tline(18:21) = '  U ';
                        cres = resnum;
                        fprintf(ofid,'%s\n',tline);
                    case 'RCYT'
                        tline(18:21) = '  C ';
                        cres = resnum;
                        fprintf(ofid,'%s\n',tline);
                    case 'ADE '
                        tline(18:21) = ' DA ';
                        cres = resnum;
                        fprintf(ofid,'%s\n',tline);
                    case 'CYT '
                        tline(18:21) = ' DC ';
                        cres = resnum;
                        fprintf(ofid,'%s\n',tline);
                    case 'GUA '
                        tline(18:21) = ' DG ';
                        cres = resnum;
                        fprintf(ofid,'%s\n',tline);
                    case 'THY '
                        tline(18:21) = ' DT ';
                        cres = resnum;
                        fprintf(ofid,'%s\n',tline);
                    case {'CYSS','CYSM'}
                        tline(18:21) = 'CYS ';
                        cres = resnum;
                        switch tline(13:16)
                            case {' C  ',' O  ',' N  ',' H  ',' CA ',' HA ',' CB ',' HB2',' HB3',' SG '}
                                fprintf(ofid,'%s\n',tline);
                            otherwise
                                atnum = atnum - 1;
                        end
                    otherwise
                        cres = resnum;
                        fprintf(ofid,'%s\n',tline);
                end
            end
        else
            is_atom_line = false;
        end
    end
    if strcmp(record,'TER   ')
        fprintf(ofid,'%s\n',tline);
    end
    if ~skipped && is_atom_line
        tline0 = tline;
    end
end

fclose(ifid);
fclose(ofid);

if ~isempty(entity)
    [entity,exceptions] = get_pdb(tmpname,options,entity);
else
    [entity,exceptions] = get_pdb(tmpname,options);
end

delete(tmpname);