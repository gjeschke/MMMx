function [entity,exceptions] = get_pdb(ident,options)
%
% GET_PDB Load structure from PDB server or local file
%
%   [entity,exceptions] = GET_PDB
%   Returns an (entity) structure in MMMx:atomic representation
%
%   [entity,exceptions] = GET_PDB(ident)
%   Reads topology and Cartesian coordinates of a biological entity
%   from a PDB file obtained from a server or locally
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
%
% OUTPUT
% entity       entity structure in MMMx:atomic representation
% exceptions   error message if something went wrong, entity is empty for
%              errors but not for warnings, for warnings, only the last one
%              is reported, cell containing an empty array, if no exception
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

profile on

% use European PDB server
PDB_server = 'ftp.ebi.ac.uk';
PDB_structures = '/pub/databases/rcsb/pdb-remediated/data/structures/all/pdb/';

% maximum number of atoms for array pre-allocation, function gets slow, if
% tghis number is too small and is memory-intensive, if it is too large
maxatoms = 1000000;
maxwater = 100000;

% initialize empty outputs
entity = [];
exceptions{1} = [];
warnings = 0; % counter for warnings

if ~exist('options','var') || ~isfield(options,'dssp') || ...
        isempty(options.dssp)
    options.dssp = false;
end

% placeholder for dowloaded file later to be deleted
to_be_deleted = ''; 


[~,~,ext] = fileparts(ident);

% load file from PDB server if required 
if isempty(ext) && length(ident) == 4
    entity.name = ident;
    fname0 = ['pdb' lower(ident) '.ent.gz'];
    % open FTP connection
    try
        ftp_obj = ftp(PDB_server,'anonymous','anonymous');
    catch exception
        exceptions{1} = exception;
        return;
    end
    % change FTP object to structure directory
    try
    cd(ftp_obj,PDB_structures);
    catch exception
        exceptions{1} = exception;
        return;
    end
    % switch FTP object to binary mode
    try
        binary(ftp_obj);
    catch exception
        exceptions{1} = exception;
        return;
    end
    % download zipped PDB structure file 
    try
        mget(ftp_obj,fname0);
    catch exception
        exceptions{1} = exception;
        return;
    end
    % unzip the downloded file
    try
        gunzip(fname0);
    catch exception
        exceptions{1} = exception;
        return;
    end
    ident = ['pdb' lower(ident) '.ent'];
    to_be_deleted = ident;
    % delete the zipped file
    try
        delete(fname0);
    catch exception
        warnings = warnings + 1;
        exceptions{warnings} = exception;
        % if everything else went well, this produces only a warning
    end
end

% open the PDB file
try
    fid = fopen(ident);
catch exception
    warnings = warnings + 1;
    exceptions{warnings} = exception;
    return
end

preserve_residue_numbers = true; % flag for preserving residue numbers 
                                 % from PDB file
                                 
% pre-allocate Cartesian coordinate, element, and index arrays
xyz = zeros(maxatoms,3);
elements = zeros(1,maxatoms,'uint8');
occupancies = zeros(1,maxatoms,'uint8');
index_array = zeros(maxatoms,5,'uint16');
water_indices = zeros(1,maxwater,'uint32');
atoms = 0;
indexed_atoms = 0;
water_atoms = 0;
residues = 0;
chains = '';
curr_chain = '';
models = 1;
current_model = 1;
while 1
    tline = fgetl(fid);
    if ~ischar(tline) 
        break 
    end
    if length(tline) >= 66 && strcmpi(tline(1:6),'HEADER')
        entity.name = tline(63:66);
    end
    if length(tline) >= 7 && strcmpi(tline(1:5),'MODEL')
        current_model = str2double(tline(7:14));
        if current_model > models
            models = current_model;
        end
    end
    % ATOM record must contain at least 54 characters (up to coordinates) 
    % atom numbers and insertion codes are ignored
    if length(tline) >= 54 && (strcmpi(tline(1:4),'ATOM') || ...
            strcmpi(tline(1:6),'HETATM')) % atom loop
        chain = tline(22);
        if ~strcmp(chain,curr_chain)
            curr_chain = chain;
            curr_resnum = 0; % current residue number
        end
        atoms = atoms + 1;
        atname = strtrim(tline(13:16));
        for kk = 1:length(atname)
            if atname(kk) == ''''
                atname(kk) = '_';
            end
        end
        altloc = tline(17);
        resname = strtrim(tline(18:20));
        if ~contains(chains,chain)
            chains = strcat(chains,chain);
        end
        if isstrprop(chain,'lower')
            chainfield = strcat(upper(chain),'_');
        else
            chainfield = chain;
        end
        trial_resnum = str2double(tline(23:26));
        if trial_resnum < curr_resnum
            preserve_residue_numbers = false;
            trial_resnum = curr_resnum + 1;
        end
        % the infamous insertion code
        if tline(27) ~= ' '
            if trial_resnum ~= curr_resnum
                trial_resnum = curr_resnum + 1;
                offset = offset + 1;
            end
            preserve_residue_numbers = false;
        end
        if trial_resnum ~= curr_resnum
            if curr_resnum > 0
                entity.(chainfield).(resfield).locations = altlocs;
            end
            altlocs = altloc;
            if altlocs ~= ' '
                altlocs = [' ' altlocs]; %#ok<AGROW>
            end
            atnames = [':' atname,':'];
        else
            if ~contains(altlocs,altloc)
                altlocs = [altlocs altloc]; %#ok<AGROW>
            end
            if ~contains(atnames,[':' atname ':'])
                atnames = strcat(atnames,strcat(atname,':'));
            end
        end 
        x = str2double(tline(31:38)); 
        y = str2double(tline(39:46)); 
        z = str2double(tline(47:54));
        % initialize all attributes that may not be present
        occupancy = 1;
        Bfactor = NaN;
        charge = 0;
        if length(tline) >= 60
            occupancy = str2double(tline(55:60));
        end
        if length(tline) >= 66
            Bfactor = str2double(tline(61:66));
        end
        if length(tline) >= 78
            element = strtrim(tline(77:78));
        else
            element = get_element_by_atomname(atname);
        end
        if length(tline) >= 80
            charge = str2double(tline(79:80));
            if isnan(charge)
                charge = 0;
            end
        end
        if trial_resnum > residues
            residues = trial_resnum;
        end
        xyz(atoms,:) = [x,y,z];
        elements(atoms) = element_number(element);
        occupancies(atoms) = round(100*occupancy);
        if strcmpi(resname,'HOH') % special treatment for water
            water_atoms = water_atoms + 1;
            water_indices(water_atoms) = atoms;
        else % atoms other than water
            resfield = sprintf('R%i',trial_resnum);
            % update topology information
            chain_index = strfind(chains,chain);
            entity.(chainfield).index = chain_index; 
            entity.(chainfield).(resfield).index = trial_resnum;
            entity.(chainfield).(resfield).name = resname;
            entity.(chainfield).(resfield).populations = 1; % single rotamer
            entity.(chainfield).(resfield).(atname).element = element;
            entity.(chainfield).(resfield).(atname).charge = charge;
            entity.(chainfield).(resfield).(atname).bfactor = Bfactor;
            % determine missing rotamer indices
            rotamer_index = strfind(altlocs,altloc);
            atom_index = strfind(atnames,[':' atname ':']);
            atom_index = 1 + sum(find(atnames == ':') < atom_index(1));
            entity.(chainfield).(resfield).(atname).index = atom_index;
            % entry into index array
            indexed_atoms = indexed_atoms + 1;
            index_array(atoms,:) = [chain_index,trial_resnum,atom_index,current_model,rotamer_index];
        end
        % update current residue number
        curr_resnum = trial_resnum;
    end % end atom loop
end
entity.original_residue_numbers = preserve_residue_numbers;
entity.xyz = xyz(1:atoms,:);
entity.elements = elements(1:atoms);
entity.occupancies = occupancies(1:atoms);
entity.index_array = index_array(1:indexed_atoms,:);
entity.water = water_indices(1:water_atoms);
% add conformer populations
entity.populations = ones(1,models)/models;

% set entity name, if it was not yet assigned
if ~isfield(entity,'name')
    entity.name = 'MMMx';
end
if isfield(options,'name') && ~isempty(options.name)
    entity.name = options.name;
end

% close the PDB file
try
    fclose(fid);
catch exception
    warnings = warnings + 1;
    exceptions{warnings} = exception;
end

if options.dssp && preserve_residue_numbers
    dssp = get_dssp(ident);
    for k = 1:length(dssp)
        if isstrprop(dssp(k).chain,'lower')
            chainfield = strcat(upper(chain),'_');
        else
            chainfield = dssp(k).chain;
        end
        resfield = strcat('R',dssp(k).tag);
        entity.(chainfield).(resfield).dssp = dssp(k).sec;
        entity.(chainfield).(resfield).sheet = dssp(k).bp;
    end
end

% delete the dowloaded file, if one was downloaded
if ~isempty(to_be_deleted)
    try
        delete(to_be_deleted);
    catch exception
        warnings = warnings + 1;
        exceptions{warnings} = exception;
    end
end

profile viewer

% functions private to get_pdb follow

function element = get_element_by_atomname(atom_tag)
% derive element from atom name, this is safe only if two-character
% elements have atom names shorter than four characters
% normally, any PDB file that has two-letter elements should also have
% specified elements in characters 77-78

element = atom_tag;
for k = 1:length(element)
    if isstrprop(element(k),'digit'), element(k)=' '; end
end
if element(1)=='H' && length(atom_tag)==4 % fix Biomer format
    element='H';
end
if element(1)=='C' && length(atom_tag)==4 % fix Biomer format
    element='C';
end
if element(1)=='N' && length(atom_tag)==4 % fix Biomer format
    element='N';
end
if element(1)=='O' && length(atom_tag)==4 % fix Biomer format
    element='O';
end
element=strtrim(element);

function elnum = element_number(element)

pse='#HHeLiBe#B#C#N#O#FNeNaMgAlSi#P#SClAr#KCaScTi#VCrMnFeCoNiCuZnGaGeAsSeBrKrRbSr#YZrNbMoTcRuRhPdAgCdInSnSbTe#IXeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTa#WReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaUNpPuAmCmBkCfEsFmMdNoLr';

if length(element) < 2
    element = strcat('#',element);
else
    element(2) = lower(element(2));
end
elnum = (strfind(pse,element)+1)/2;

function dssp = get_dssp(fname)

dssp=[];

dospath = which('dssp-2.0.4-win32.exe');
if isempty(dospath)
    dospath = which('dssp.exe');
end
if isempty(dospath)
    dospath = which('dsspcmbi.exe');
end
if ~isempty(dospath) % suppress this if DSSP not known
    infile = which(fname);
    poi = strfind(infile,'.pdb');
    if isempty(poi)
        poi = strfind(infile,'.ent');
    end
    if isempty(poi) || poi(1)<2
        basname = infile;
    else
        basname = infile(1:poi(1)-1);
    end
    dssp_file = [basname '.dssp'];
    cmd=[dospath ' ' infile ' ' dssp_file];
    s = dos(cmd);
    if s~=0 % dssp was not successful
        return
    else
        dssp = rd_dssp(dssp_file);
    end
end 


