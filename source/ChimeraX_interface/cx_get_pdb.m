function [id,entity] = cx_get_pdb(ident,options)
%
% CX_GET_PDB Load PDB structure into ChimeraX and retrieve entity
%
%   id = CX_GET_PDB(ident)
%   Loads entity into ChimeraX and provides the ChimeraX identifier
%
%   [id,entity] = CX_GET_PDB(ident)
%   Returns ChimeraX model number and entity in MMMx:atomic representation
%
%   [id,entity] = CX_GET_PDB(ident,options)
%   Loads structure from PDB server or local file into ChimeraX vi its
%   'open' command and returns the structure as an entity in MMMx:atomic
%   format, for future reference, the ChimeraX identifier is returned
%
% INPUT
% ident     PDB identifier, if four characters without extension, otherwise
%           file name, can be 'browser' or empty, in that case a ChimeraX
%           browser window is opened
% options   structure with requests for extended information, flags
%           .dssp    DSSP (Kabsch/Sanders) secondary structure
%                    defaults to false
%           .name    optional name for the entity, defaults: PDB identifier
%                    upon download or ChimeraX name
%           .altlocs alternate locations
%           .Bfactor B factors
%           .element retrieve from ChimeraX, otherwise derived from atom
%                    name
%           .MMMcol  apply MMM color scheme to ribbon, default true
%
% OUTPUT
% id           internal identifier of ChimeraX
% entity       entity structure in MMMx:atomic representation
%
% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke


% maximum number of atoms for array pre-allocation, function gets slow, if
% this number is too small and is memory-intensive, if it is too large
maxatoms = 1000000;
maxwater = 100000;

id = [];
entity = [];

if ~exist('ident','var') || isempty(ident)
    ident = 'browse';
end

if ~exist('options','var') || ~isfield(options,'dssp') || ...
        isempty(options.dssp)
    options.dssp = false;
end

if ~isfield(options,'altlocs') || isempty(options.altlocs)
    options.altlocs = false;
end

if ~isfield(options,'Bfactor') || isempty(options.Bfactor)
    options.Bfactor = false;
end

if ~isfield(options,'element') || isempty(options.element)
    options.element = false;
end

if ~isfield(options,'MMMcol') || isempty(options.MMMcol)
    options.MMMcol = true;
end

[path,basname,ext] = fileparts(ident);

if ~isempty(path)
    response = cx_command(sprintf('cd %s',path));
    fprintf(1,'%s\n',response);
elseif ~isempty(ext) || (length(ident) ~= 4 && ~strcmp(ident,'browse'))
    response = cx_command(sprintf('cd %s',pwd));
    fprintf(1,'%s\n',response);
end


response = cx_command(sprintf('open %s%s',basname,ext));

poi = strfind(response,'Chain information for');
if isempty(poi)
    return
end
idstr= split(response(poi:end),'#');
idstr = strtok(idstr);
idstr = idstr{2};
poi = strfind(idstr,'.');
if ~isempty(poi)
    idstr = idstr(1:poi-1);
end
id = str2double(idstr);
if isempty(id) || isnan(id)
    return
end

if options.MMMcol
    cx_command(sprintf('color #%i %s',id,'light goldenrod yellow target c'));
    cx_command(sprintf('color #%i %s',id,'& helix lightsalmon target c'));
    cx_command(sprintf('color #%i %s',id,'& strand steelblue target c'));
end

if nargout > 1
    my_info = cx_get_attribute(sprintf('#%i',id),'model','name');

    if contains(lower(my_info),'group"')
        poi = [strfind(my_info,'model id') length(my_info)];
        model = 0;
        model_list = cell(1,length(poi)-1);
        for mod = 1:length(poi)-1
            modstr = my_info(poi(mod):poi(mod+1)-2);
            if contains(modstr,'AtomicStructure')
                args = split(modstr);
                model = model+1;
                model_list{model} = args{3};
                my_name = args{7};
            end
        end
        model_list = model_list(1:model);
    else
        my_name = my_info(1:4);
        model_list{1} = sprintf('#%i',id);
    end
    if ~isfield(options,'name') || isempty(options.name)
        entity.name = upper(my_name);
    else
        entity.name = options.name;
    end
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
    models = length(model_list);
    current_model = 1;
    old_resname = 'HOH'; % avoid location entry for water residues
    populations = zeros(length(models),1);
    conformers = 0;
    for mod = 1:length(model_list)
        curr_chain = '';
        response = cx_command(sprintf('getcrd %s',model_list{mod}));
        atoms_info = splitlines(response);
        for k = 1:length(atoms_info)
            if ~isempty(atoms_info{k})
                clear locs
                curr_atom = split(atoms_info{k});
                address = cx_to_mmmx_address(curr_atom{2});
                altlocs = ' ';
                if ~strcmp(address.chain,curr_chain)
                    curr_chain = address.chain;
                    curr_resnum = 0; % current residue number
                end
                atname = strtrim(upper(address.atom));
                for kk = 1:length(atname)
                    if atname(kk) == ''''
                        atname(kk) = '_';
                    end
                end
                if ~isstrprop(atname(1),'alpha')
                    atname = strcat('Z_',atname);
                end
                resnum_str = strtrim(address.residue);
                if ~contains(chains,address.chain)
                    chains = strcat(chains,address.chain);
                end
                if isstrprop(address.chain,'lower')
                    chainfield = strcat(upper(address.chain),'_');
                else
                    chainfield = address.chain;
                end
                trial_resnum = str2double(resnum_str);
                if trial_resnum ~= curr_resnum
                    resname_info = cx_get_attribute(curr_atom{2},'residue','name');
                    resname_str = split(resname_info);
                    resname = resname_str{1};
                end
                if options.altlocs
                    % check whether this atom has alternate locations
                    loc_info = cx_get_attribute(curr_atom{2},'atom','alt_locs');
                    if ~isempty(loc_info)
                        locs = split(loc_info,',');
                        for l = 1:length(locs)
                            if ~contains(altlocs,locs{l})
                                altlocs = [altlocs,locs{l}]; %#ok<AGROW>
                            end
                        end
                    else
                        locs{1} = '';
                    end
                else
                    locs{1} = '';
                end
                if trial_resnum ~= curr_resnum
                    if curr_resnum > 0 && ~strcmpi(old_resname,'HOH')
                        entity.(chainfield).(resfield).locations = altlocs;
                    end
                    atnames = [':' atname,':'];
                else
                    if ~contains(atnames,[':' atname ':'])
                        atnames = strcat(atnames,strcat(atname,':'));
                    end
                end
                occupancy = 1;
                Bfactor = 0;
                charge = 0;
                for loc = 1:length(locs)
                    atoms = atoms + 1;
                    if isempty(locs{loc})
                        x = str2double(curr_atom{3});
                        y = str2double(curr_atom{4});
                        z = str2double(curr_atom{5});
                    else
                        cx_set_attribute(curr_atom{2},'atom','alt_loc',locs{loc});
                        coor_info = cx_get_attribute(curr_atom{2},'atom','coord');
                        coor_str = split(coor_info,',');
                        x = str2double(coor_str{1});
                        y = str2double(coor_str{2});
                        z = str2double(coor_str{3});
                    end
                    if options.Bfactor
                        Bfactor_str = cx_get_attribute(curr_atom{2},'atom','bfactor');
                        Bfactor =str2double(Bfactor_str);
                    end
                    if options.element
                        element = cx_get_attribute(curr_atom{2},'atom','element');
                    else
                        element = get_element_by_atomname(atname);
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
                    else
                        resfield = sprintf('R%i',trial_resnum);
                        % update topology information
                        chain_index = strfind(chains,address.chain);
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
                        rotamer_index = strfind(altlocs,locs{loc});
                        if isempty(rotamer_index)
                            rotamer_index = 1;
                        end
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
                end
            end
        end
    end
    populations = populations(1:conformers);
    entity.original_residue_numbers = true;
    entity.xyz = xyz(1:atoms,:);
    entity.elements = elements(1:atoms);
    entity.occupancies = occupancies(1:atoms);
    entity.index_array = index_array(1:indexed_atoms,:);
    entity.water = water_indices(1:water_atoms);
    entity.water_selected = false;
    % add conformer populations
    if conformers == models
        entity.populations = populations;
    else
        entity.populations = ones(models,1)/models;
    end
    entity.selected = 1; % first conformer is selected by default
end

% functions private to cx_get_pdb follow

function element = get_element_by_atomname(atom_tag)
% derive element from atom name, this is safe only if two-character
% elements have atom names shorter than four characters
% normally, any PDB file that has two-letter elements should also have
% specified elements in characters 77-78

element = atom_tag;
for k = 1:length(element)
    if isstrprop(element(k),'digit'), element(k)=' '; end
end
if element(1)=='H' 
    element='H';
end
if element(1)=='C' 
    element='C';
end
if element(1)=='N' 
    element='N';
end
if element(1)=='O' 
    element='O';
end
if element(1)=='S' 
    element='S';
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

