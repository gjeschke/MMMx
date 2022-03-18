function [entity,exceptions] = get_pdb(ident,options,entity)
%
% GET_PDB Load structure from PDB server or local file
%
%   [entity,exceptions] = GET_PDB(ident)
%   Returns an (entity) structure in MMMx:atomic representation
%
%   [entity,exceptions] = GET_PDB(ident,options)
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
%           .stripH Boolean flag indcating that protons should be removed,
%                   that are not expected to be always attached
%                   needed to generate consistent ensembles, defaults to
%                   false
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
% Copyright(c) 2020: Gunnar Jeschke

% maximum number of atoms for array pre-allocation, function gets slow, if
% this number is too small and is memory-intensive, if it is too large
maxatoms = 1000000;
maxwater = 100000;
maxmodels = 100000;
resdefs = load('monomers.mat');

% initialize empty outputs
if ~exist('entity','var')
    entity = [];
    add_to_entity = false;
else
    add_to_entity = true;
end
exceptions{1} = [];
warnings = 0; % counter for warnings

if ~exist('options','var') || isempty(options) || ~isfield(options,'dssp') || ...
        isempty(options.dssp)
    options.dssp = false;
end
if ~isfield(options,'stripH') || isempty(options.stripH)
    options.stripH = false;
elseif options.stripH
    min_atoms = load('minimal_protonation.mat');
end

% placeholder for downloaded file later to be deleted
to_be_deleted = ''; 


[~,~,ext] = fileparts(ident);

% load file from PDB server if required
if isempty(ext) && length(ident) == 4
    query = sprintf('https://files.rcsb.org/download/%s.pdb.gz',lower(ident));
    fname0 = [lower(ident) '.pdb.gz'];
    try
        websave(fname0,query);
    catch exception
        exceptions{1} = exception;
        return;
    end
    try
        gunzip(fname0);
    catch exception
        exceptions{1} = exception;
        return;
    end
    ident = sprintf('%s.pdb',lower(ident));
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
occupancies = zeros(maxatoms,1,'uint8');
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
model_offset = 0;
populations = ones(maxmodels,1);
conformers = 0;
if ~isempty(entity)
    model_offset = length(entity.populations);
    current_model = model_offset + 1; 
    atoms = length(entity.elements);
    elements(1:atoms) = entity.elements;
    occupancies(1:atoms) = entity.occupancies;
    [indexed_atoms,~] = size(entity.index_array);
    index_array(1:indexed_atoms,:) = entity.index_array;
    xyz(1:atoms,:) = entity.xyz;
    water_atoms = length(entity.water);
    water_indices(1:water_atoms) = entity.water;
    conformers = length(entity.populations) + 1;
    populations(1:conformers-1) = entity.populations;
    populations(conformers) = mean(populations);
end
old_resname = 'HOH'; % avoid location entry for water residues
% initialize information on modified residues
modres(1000).actual = '';
modres(1000).clear_name = '';
modres(1000).native = '';
modres(1000).chain = '';
modres(1000).residue = [];
modified_residues = 0;
offset = 0; % residue number offset to avoid negative residue numbers

while 1
    tline = fgetl(fid);
    if ~ischar(tline) 
        break 
    end
    if length(tline) >= 66 && strcmpi(tline(1:6),'HEADER')
        entity.name = tline(63:66);
    end
    if length(tline) >= 30 && strcmpi(tline(1:6),'MODRES')
        modified_residues = modified_residues + 1;
        modres(modified_residues).actual = tline(13:15);
        modres(modified_residues).clear_name = strtrim(tline(30:end));
        modres(modified_residues).native = tline(25:27); 
        modres(modified_residues).chain = tline(17); 
        modres(modified_residues).residue = str2double(tline(18:22)); 
    end
    if length(tline) >= 7 && strcmpi(tline(1:5),'MODEL')
        if length(tline) < 14
            tline = pad(tline,14); % because PED does not check for correct format
        end
        current_model = str2double(tline(7:14)) + model_offset;
        if current_model > models
            models = current_model;
        end
        curr_chain = '';
    end
    % read population information in MMMx:atomic PDB files
    if length(tline) >= 48 && contains(tline,'REMARK 400   MODEL') && contains(tline,'POPULATION')
        conformer = str2double(tline(19:28)) + model_offset;
        if conformer > 0
            populations(conformer) = str2double(tline(40:48));
            if conformer > conformers
                conformers = conformer;
            end
        end
    end   
    % ATOM record must contain at least 54 characters (up to coordinates) 
    % atom numbers and insertion codes are ignored
    if length(tline) >= 54 && (strcmpi(tline(1:4),'ATOM') || ...
            strcmpi(tline(1:6),'HETATM')) % atom loop
        chain = tline(22);
        if chain == ' '
            chain = 'Z';
        end
        if ~strcmp(chain,curr_chain)
            curr_chain = chain;
            curr_resnum = 0; % current residue number
            offset = 0; % residue number offset to avoid negative residue numbers
        end
        atname0 = tline(13:16);
        atname = strtrim(upper(tline(13:16)));
        for kk = 1:length(atname)
            if atname(kk) == ''''
                atname(kk) = '_';
            end
            if atname(kk) == '"'
                atname(kk) = '_';
            end
        end
        if ~isstrprop(atname(1),'alpha')
            atname = strcat('Z_',atname);
        end
        if length(tline) >= 78
            element = strtrim(tline(77:78));
        else
            element = get_element_by_atomname(atname0);
        end
        elm_num = element_number(element);
        resname = strtrim(tline(18:20));
        if options.stripH && elm_num == 1
            canonical = find_proton(resname,atname,min_atoms);
            if ~canonical
                continue
            end
        end
        atoms = atoms + 1;
        altloc = tline(17);
        if ~contains(chains,chain)
            chains = strcat(chains,chain);
        end
        if isstrprop(chain,'lower')
            chainfield = strcat(upper(chain),'_');
        else
            chainfield = chain;
        end
        trial_resnum = str2double(tline(23:26));
        if trial_resnum < 1 && offset ==0
            offset = 1 - trial_resnum;
        end
        trial_resnum = trial_resnum + offset;
%         if trial_resnum < curr_resnum
%             preserve_residue_numbers = false;
%             trial_resnum = curr_resnum + 1;
%         end
        % the infamous insertion code
        if tline(27) ~= ' '
            if trial_resnum ~= curr_resnum
                trial_resnum = curr_resnum + 1;
                offset = offset + 1;
            end
            preserve_residue_numbers = false;
        end
        if trial_resnum ~= curr_resnum
            if curr_resnum > 0 && ~strcmpi(old_resname,'HOH')
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
        [elements(atoms),element] = element_number(element);
        occupancies(atoms) = round(100*occupancy);
        if strcmpi(resname,'HOH') % special treatment for water
            water_atoms = water_atoms + 1;
            water_indices(water_atoms) = atoms;
        else % atoms other than water
            resfield = sprintf('R%i',trial_resnum);
            % update topology information
            chain_index = strfind(chains,chain);
            entity.(chainfield).index = chain_index;
            entity.(chainfield).selected = false;
            entity.(chainfield).(resfield).index = trial_resnum;
            entity.(chainfield).(resfield).selected = 0; 
            % first rotamer is selected by default
            entity.(chainfield).(resfield).selected_rotamers = 1;
            entity.(chainfield).(resfield).name = resname;
            entity.(chainfield).(resfield).locations = altlocs;
            entity.(chainfield).(resfield).populations = 1; % single rotamer
            entity.(chainfield).(resfield).(atname).element = element;
            entity.(chainfield).(resfield).(atname).charge = charge;
            entity.(chainfield).(resfield).(atname).bfactor = Bfactor;
            entity.(chainfield).(resfield).(atname).selected = 0;
            entity.(chainfield).(resfield).(atname).selected_locations = 1;
            if ~isfield(entity.(chainfield).(resfield).(atname),'tab_indices')
                entity.(chainfield).(resfield).(atname).tab_indices = atoms;
            else
                entity.(chainfield).(resfield).(atname).tab_indices = ...
                    [entity.(chainfield).(resfield).(atname).tab_indices atoms];
            end
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
        old_resname = resname;
%         if atoms == 992
%             disp('Aber hallo!');
%         end
    end % end atom loop
end
if conformers == 0
    conformers = 1;
    populations(1) = 1;
end
populations = populations(1:conformers);
entity.original_residue_numbers = preserve_residue_numbers;
entity.xyz = xyz(1:atoms,:);
entity.elements = elements(1:atoms);
entity.occupancies = occupancies(1:atoms);
entity.index_array = index_array(1:indexed_atoms,:);
entity.water = water_indices(1:water_atoms);
entity.water_selected = false;
% add conformer populations
if conformers == models || add_to_entity
    entity.populations = populations;
else
    entity.populations = ones(models,1)/models;
end
entity.selected = 1; % first conformer is selected by default

% set entity name, if it was not yet assigned
if ~isfield(entity,'name')
    entity.name = 'MMMx';
end
if isfield(options,'name') && ~isempty(options.name)
    entity.name = options.name;
end

% add information on modified resisues
modres = modres(1:modified_residues);
entity.modres = modres;


% close the PDB file
try
    fclose(fid);
catch exception
    warnings = warnings + 1;
    exceptions{warnings} = exception;
end

if options.dssp && preserve_residue_numbers
    dssp = get_dssp(ident);
    % the following is required, since DSSP sometimes defines new chains
    chainfields = fieldnames(entity);
    cpoi = 0;
    for kc = 1:length(chainfields)
        chain = chainfields{kc};
        if isstrprop(chain(1),'upper')
            cpoi = cpoi + 1;
            chainfields{cpoi} = chainfields{kc};
        end
    end
    chainfields = chainfields(1:cpoi);
    for k = 1:length(dssp)
        if isstrprop(dssp(k).chain,'lower')
            chainfield = strcat(upper(chain),'_');
        else
            chainfield = dssp(k).chain;
        end
        chain_exists = false;
        for kc = 1:length(chainfields)
            if strcmp(chainfield,chainfields{kc})
                chain_exists = true;
            end
        end
        if chain_exists
            resfield = strcat('R',dssp(k).tag);
            entity.(chainfield).(resfield).dssp = dssp(k).sec;
            entity.(chainfield).(resfield).sheet = dssp(k).bp;
        end
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
elseif element(1)==' ' && element(2) == 'H' % fix CYANA
    element = 'H';
end
if element(1)=='C' && length(atom_tag)==4 % fix Biomer format
    element='C';
elseif element(1)==' ' && element(2) == 'C'
    element = 'C'; % fix CYANA
end
if element(1)=='N' && length(atom_tag)==4 % fix Biomer format
    element='N';
elseif element(1)==' ' && element(2) == 'N'
    element = 'N'; % fix CYANA
end
if element(1)=='O' && length(atom_tag)==4 % fix Biomer format
    element='O';
elseif element(1)==' ' && element(2) == 'O'
    element = 'O'; % fix CYANA
end
element=strtrim(element);

function [elnum,element] = element_number(element)

pse='#HHeLiBe#B#C#N#O#FNeNaMgAlSi#P#SClAr#KCaScTi#VCrMnFeCoNiCuZnGaGeAsSeBrKrRbSr#YZrNbMoTcRuRhPdAgCdInSnSbTe#IXeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTa#WReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaUNpPuAmCmBkCfEsFmMdNoLr';

if length(element) < 2
    element = strcat('#',element);
else
    element(2) = lower(element(2));
end
elnum = (strfind(pse,element)+1)/2;
if isempty(elnum)
    element = element(1);
    elnum = (strfind(pse,['#',element])+1)/2;
end
if element(1) == '#'
    element = element(2:end);
end

function canonical = find_proton(resname,atname,min_atoms)
% checks whether a proton is canonical, if the residue is not an amino
% acid or nucleic acid, the proton is assumed to be canonical

canonical = 0;
if isfield(min_atoms.minimal_atoms,resname)
    pattern = fieldnames(min_atoms.minimal_atoms.(resname));
    for kp = 1:length(pattern)
        this_atom = pattern{kp};
        if strcmp(atname,this_atom)
            canonical = min_atoms.minimal_atoms.(resname).(this_atom);
        end
    end
else
    canonical = 1;
end