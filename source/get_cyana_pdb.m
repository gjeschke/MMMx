function [entity,exceptions] = get_cyana_pdb(fname,options,entity)
%
% GET_CYANA_PDB Load structure from CYANA pseudo-PDB file
%
%   [entity,exceptions] = GET_CYANA_PDB(fname)
%   Returns an (entity) structure in MMMx:atomic representation
%
%   [entity,exceptions] = GET_CYANA_PDB(fname,options)
%   Reads topology and Cartesian coordinates of a biological entity
%   from a CYANA pseudo-PDB file
%
%   [entity,exceptions] = GET_CYANA_PDB(fname,options,entity)
%   Adds Cartesian coordinates of a conformer to an input biological entity
%   by reading from a CYANA pseudo-PDB file
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

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
%
% G. Jeschke, 2021

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

skipmode = true;

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
plnum = 0;
ll5num = 0;
lnnum = 0;
orinum = 0;
nlnum = 0;
lpnum = 0;
while 1 && chain <= maxchain
    tline = fgetl(ifid);
    if ~ischar(tline), break, end
    if length(tline)<26 % catches too short end line of MolProbity files
        fprintf(ofid,'%s\n',tline);
    else
        resnum = str2double(tline(23:26));
        if resnum > 531
            tline(22) = 'B';
        end
        record = tline(1:6);
        if strcmp(record,'ATOM  ') && skipmode
            chainid = tline(22);
            resnum = str2double(tline(23:26));
            skip = true;
            if chainid == 'A' && resnum >= 58 && resnum <= 531
                skip = false;
            end
            if chainid == 'B' && resnum >= 288+315 && resnum <= 370+315
                skip = false;
            end
            if skip
                continue;
            end
        end
        restype = tline(18:21);
        switch restype
            case ' PL '
                skipped = true;
            case 'LL5 '
                skipped = true;
            case ' LN '
                skipped = true;
            case ' NL '
                skipped = true;
            case ' LP '
                skipped = true;
            case 'ORI '
                skipped = true;
            otherwise
                skipped = false;
        end
        % fprintf(1,'%s\n',record);
        if strcmp(record,'ATOM  ') || strcmp(record,'HETATM')
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
            resnum0 = str2double(tline(24:26));
            if tline(22) == 'B'
                tline(24:26) = sprintf('%i',resnum0-315);
            end
            if strcmp(record,'HETATM')
                cres = resnum;
                fprintf(ofid,'%s\n',tline);
            else % wrong ATOM records still need to be fixed
                restype = tline(18:21);
                switch restype
                    case ' PL '
                        plnum = plnum + 1;
                        atnum = atnum - 1;
                    case 'LL5 '
                        ll5num = ll5num + 1;
                        atnum = atnum - 1;
                    case ' LN '
                        lnnum = lnnum + 1;
                        atnum = atnum -1;
                    case ' NL '
                        nlnum = nlnum + 1;
                        atnum = atnum -1;
                    case ' LP '
                        lpnum = lpnum + 1;
                        atnum = atnum -1;
                    case 'ORI '
                        orinum = orinum + 1;
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
                    otherwise
                        cres = resnum;
                        fprintf(ofid,'%s\n',tline);
                end
            end
        end
    end
    if ~skipped
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