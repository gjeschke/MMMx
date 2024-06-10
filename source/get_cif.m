function entity = get_cif(fname)
%
% GET_CIF Generate entity from a local mmcif file
%
%   entity = GET_CIF(fname)
%   Returns an (entity) structure in MMMx:atomic representation
%
% INPUT
% fname     name of a local MMCIF file, extension .cif is appended if none
%
% OUTPUT
% entity       entity structure in MMMx:atomic representation
%
% REMARKS
% - uses rd_mmcif for extracting the 3D structure information from an MMCIF
%   file
% - for more complicated entity building, recycle to a PDB file and use
%   get_pdb

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2024: Gunnar Jeschke

% maximum number of atoms for array pre-allocation, function gets slow, if
% this number is too small and is memory-intensive, if it is too large
maxwater = 1000000;

% initialize empty output
entity = [];

[~,~,ext] = fileparts(fname);

% append extension if none
if isempty(ext) 
    fname = strcat(fname,'.cif');
end

% read 3D structure information from file
info = rd_mmcif(fname);
atoms = length(info.atomsite);
                                 
% pre-allocate Cartesian coordinate, element, and index arrays
xyz = zeros(atoms,3);
elements = zeros(1,atoms,'uint8');
occupancies = zeros(atoms,1,'uint8');
index_array = zeros(atoms,5,'uint16');
water_indices = zeros(1,maxwater,'uint32');
indexed_atoms = 0;
water_atoms = 0;
chains = '';
conformers = 0;
current_model = 1;
curr_resnum = 0;
old_resname = 'HOH'; % avoid location entry for water residues
clear_names(10000).chain = '';
clear_names(10000).residue = '';
clear_names(10000).atom = '';
clear_names(10000).cid = 0;
clear_names(10000).rid = 0;
clear_names(10000).aid = 0;
clear_names(10000).rotid = 0;

old_model = 0;
for a = 1:length(info.atomsite)
    x = str2double(info.atomsite(a).Cartnx);
    y = str2double(info.atomsite(a).Cartny);
    z = str2double(info.atomsite(a).Cartnz);
    xyz(a,:) = [x,y,z];
    this_model = str2double(info.atomsite(a).pdbxPDBmodelnum);
    if this_model > conformers
        conformers = this_model;
    end
    if this_model ~=  old_model
        this_atom = 0;
        old_model = this_model;
    end
    if this_model == 1
        chain = info.atomsite(a).labelasymid;
        if chain == ' '
            chain = 'Z';
        end
        if ~contains(chains,chain)
            chains = strcat(chains,chain);
            new_chain = true;
        else
            new_chain = false;
        end
        atname = info.atomsite(a).labelatomid;
        atname = strip(atname,'"');
        for kk = 1:length(atname)
            if atname(kk) == ''''
                atname(kk) = '_';
            end
            if atname(kk) == '"'
                atname(kk) = '_';
            end
            if atname(kk) == '*'
                atname(kk) = '_';
            end
        end
        if ~isstrprop(atname(1),'alpha')
            atname = strcat('Z_',atname);
        end
        element = info.atomsite(a).typesymbol;
        resname = info.atomsite(a).labelcompid;
        altloc = info.atomsite(a).labelaltid;
        if altloc == '.'
            altloc = ' ';
        end
        if isstrprop(chain,'lower')
            chainfield = strcat(upper(chain),'_');
        else
            chainfield = chain;
        end
        resnum = str2double(info.atomsite(a).labelseqid);
        if resnum ~= curr_resnum && new_chain
            altlocs = altloc;
        end
        if resnum ~= curr_resnum 
            atnames = [':' atname,':'];
            if ~new_chain
                if curr_resnum > 0 && ~strcmpi(old_resname,'HOH')
                    entity.(chainfield).(resfield).locations = altlocs;
                end
                altlocs = altloc;
            end
        else
            if ~contains(altlocs,altloc)
                altlocs = [altlocs altloc]; %#ok<AGROW>
            end
            if ~contains(atnames,[':' atname ':'])
                atnames = strcat(atnames,strcat(atname,':'));
            end
        end
        % initialize all attributes that may not be present
        occupancy = str2double(info.atomsite(a).occupancy);
        Bfactor = str2double(info.atomsite(a).Bisoorequiv);
        if isfield(info.atomsite(a),'pdbxformalcharge')
            charge = str2double(info.atomsite(a).pdbxformalcharge);
        else
            charge = 0;
        end
        if isnan(charge)
            charge = 0;
        end
        if isempty(element_number(element))
            element = 'C';
        end
        [elements(a),element] = element_number(element);
        occupancies(a) = round(100*occupancy);
        if strcmpi(resname,'HOH') % special treatment for water
            water_atoms = water_atoms + 1;
            water_indices(water_atoms) = a;
        else % atoms other than water
            this_atom = this_atom + 1;
            resfield = sprintf('R%i',resnum);
            % update topology information
            chain_index = strfind(chains,chain);
            entity.(chainfield).index = chain_index;
            entity.(chainfield).selected = false;
            entity.(chainfield).(resfield).index = resnum;
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
            clear_names(this_atom).chain = chainfield;
            clear_names(this_atom).residue = resfield;
            clear_names(this_atom).atom = atname;
            if ~isfield(entity.(chainfield).(resfield).(atname),'tab_indices')
                entity.(chainfield).(resfield).(atname).tab_indices = a;
            else
                entity.(chainfield).(resfield).(atname).tab_indices = ...
                    [entity.(chainfield).(resfield).(atname).tab_indices a];
            end
            % determine missing rotamer indices
            rotamer_index = strfind(altlocs,altloc);
            atom_index = strfind(atnames,[':' atname ':']);
            atom_index = 1 + sum(find(atnames == ':') < atom_index(1));
            entity.(chainfield).(resfield).(atname).index = atom_index;
            % entry into index array
            indexed_atoms = indexed_atoms + 1;
            index_array(a,:) = [chain_index,resnum,atom_index,current_model,rotamer_index];
            clear_names(this_atom).cid = chain_index;
            clear_names(this_atom).rid = resnum;
            clear_names(this_atom).aid = atom_index;
            clear_names(this_atom).rotid = rotamer_index;
            clear_names(this_atom).element = elements(a);
        end
        % update current residue number
        curr_resnum = resnum;
        old_resname = resname;
    else % not the first model
        if strcmpi(resname,'HOH') % special treatment for water
            water_atoms = water_atoms + 1;
            water_indices(water_atoms) = a;
        else
            this_atom = this_atom + 1;
        end
        name = clear_names(this_atom);
        entity.(name.chain).(name.residue).(name.atom).tab_indices = ...
            [entity.(name.chain).(name.residue).(name.atom).tab_indices a];
        indexed_atoms = indexed_atoms + 1;
        index_array(a,:) = [name.cid,name.rid,name.aid,current_model,name.rotid];
        elements(a) = name.element;
    end
end
entity.original_residue_numbers = true;
entity.xyz = xyz(1:atoms,:);
entity.elements = elements(1:atoms);
entity.occupancies = occupancies(1:atoms);
entity.index_array = index_array(1:indexed_atoms,:);
entity.water = water_indices(1:water_atoms);
entity.water_selected = false;
entity.clear_names = clear_names(1:this_atom);
% add conformer populations
entity.populations = ones(conformers,1)/conformers;
entity.selected = 1; % first conformer is selected by default
entity.name = info.entry_id;


% functions private to get_cif follow

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

