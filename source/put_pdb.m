function [exceptions,fid] = put_pdb(entity,fname,options)
%
% PUT_PDB Writes a minimal PDB file (one header line and coordinates)
%
%   exceptions = PUT_PDB(entity)
%   Writes minimal PDB file of the whole entity with the entity name
%   defining the file name
%
%   exceptions = GET_PDB(entity,fname)
%   Writes minimal PDB file of the entity with user-defined file name
%
%   exceptions = GET_PDB(entity,fname,options)
%   Writes minimal PDB file with additional specifications
%
%   [exceptions,fid] = put_pdb(entity,fname,options)
%   Does not end and does not close the file and returns file identifier
%
% INPUT
% entity    entity structure in MMMx:atomic representation
% fname     file name, '.pdb' is appended if it is missing
% options   structure controlling what is output, all fields are optional
%           .selected  if true, only selected objects are output, defaults
%                      to false
%           .chainIDs  renames chains in PDB file as compared to entity,
%                      cell(C,2) for C chains to be renamed, cell entries
%                      are strings
%                      .chain_IDs{c,1}   identifier of chain c in entity
%                      .chain_IDs{c,2}   identifier of chain c in PDB file,
%                                        must be a unique single letter
%           .charged   if true, atom charge is output, defaults to false
%           .bfactor   if true, B factor is output, defaults to false 
%                      (write zero B factor)
%           .order     allows to reorder M conformers, (1,M) double, only
%                      the M specified conformers are written in this order
%                      this overrides conformer selection in the entity,
%                      even if .selected is present and true
%           .pop       flag specifying whether conformer populations are
%                      written in REMARK 400, defaults to false
%           .pdbid     (pseudo-)PDB identifier, four characters, defaults 
%                      to the first four characters of entity.name, 
%                      if entity.name is empty or shorter, it is MMMX
%                      the orginal entity name is written as REMARK 400
%           .hetframe  flag, if true, write out only HETATM records for
%                      SCWRL4 heteroatom frame, defaults to false
%           .cleaned   flag, if true, write out residues without side 
%                      chains for SCWRL4 repacking, defaults to false
%           .dxyz      coordinate offset, double(1,3), in Angstroem,
%                      defaults to [0,0,0]]
%           .datnum    atom number offset, defaults to zero
%           .dresnum   residue number offset, defaults to zero
%           .fid       if present and not -1, output is appended to the
%                      file with identifier fid
%           
%
% OUTPUT
% exceptions   error message if something went wrong, cell containing an 
%              empty array, if no exception
% fid          if requested, the END record is not written, the file is not
%              closed, and the file identifier is returned
%
% the PDB identifier in the HEADER line is derived from the first four
% characters of entity.name, if entity.name is empty or shorter, it is MMMX
% the orginal entity name
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty outputs
exceptions{1} = [];

% load the monomer definitions for atom typing
biomer = load('monomer_attributes.mat');

% initialize missing input
if ~exist('fname','var') || isempty(fname)
    fname = entity.name;
    if isempty(fname)
        fname = 'MMMx';
    end
end

if ~exist('options','var') || ~isfield(options,'selected') || ...
        isempty(options.selected)
    options.selected = false;
end

if ~isfield(options,'charge') || isempty(options.charge)
    options.charge = false;
end

if ~isfield(options,'hetframe') || isempty(options.hetframe)
    options.hetframe = false;
end

if ~isfield(options,'cleaned') || isempty(options.cleaned)
    options.cleaned = false;
end

if ~isfield(options,'bfactor') || isempty(options.bfactor)
    options.bfactor = false;
end

if ~isfield(options,'dxyz') || isempty(options.dxyz)
    options.dxyz = zeros(1,3);
end

if ~isfield(options,'datnum') || isempty(options.datnum)
    options.datnum = 0;
end

if ~isfield(options,'dresnum') || isempty(options.dresnum)
    options.dresnum = 0;
end

if ~isfield(options,'fid') || isempty(options.fid)
    options.fid = -1;
end

if ~isfield(options,'order')
    options.order = [];
end

% apply coordinate offset
[na,~] = size(entity.xyz);
entity.xyz = entity.xyz + repmat(options.dxyz,na,1);

% complete file name, if neceesary
if ~contains(fname,'.pdb')
    fname = strcat(fname,'.pdb');
end

% generate (pseudo-)PDB ID, if not provided
if ~isfield(options,'pdbid') || isempty(options.pdbid) || length(options.pdbid) ~= 4
    if length(entity.name) >= 4
        pdbid = entity.name(1:4);
    else
        pdbid = 'MMMX';
    end
else
    pdbid = options.pdbid;
end

% check and determine conformer order
all_conformers = false;
if isfield(options,'order') && ~isempty(options.order)
    conformer_order = options.order;
    all_conformers = true;
    if max(conformer_order) > length(entity.populations)
        exceptions = {MException('put_pdb:unknown_conformer', 'Request to write unknown conformer %i',max(conformer_order))};
        return
    end
else
    conformer_order = 1:length(entity.populations);
end

% check for chain identifier replacement request
if isfield(options,'chainIDs') && ~isempty(options.chainIDs)
    [replace_chains,~] = size(options.chainIDs);
else
    replace_chains = 0;
end

if options.fid == -1
    [fid,message] = fopen(fname,'wt');
else % append to existing file
    fid = options.fid;
end

if fid == -1 % file could not be opened
    exceptions = {MException('put_pdb:file_could_not_be_opened', 'File could not be opened (%s)',message)};
    return
end

if options.fid == -1
    % write header line
    fprintf(fid,'HEADER    MMMX MINIMAL PDB FILE                   %s   %s\n', ...
        datestr(datetime,'dd-mmm-yy'),pdbid);

    % write original entity name as remark, lines are padded to 80 characters
    pdbline = 'REMARK 400 COMPOUND';
    fprintf(fid,'%s\n',pad(pdbline,80));
    pdbline = sprintf('REMARK 400 ENTITY NAME: %s',entity.name);
    % the line should not have more than 80 characters
    if length(pdbline) > 80
        pdbline = pdbline(1:80);
    end
    fprintf(fid,'%s\n',pad(pdbline,80));
end

% write conformer population information if requested
if isfield(options,'pop') && ~isempty(options.pop) && options.pop
    pdbline = 'REMARK 400  POPULATIONS';
    fprintf(fid,'%s\n',pad(pdbline,80));
    save_conf = 0;
    for kconf = 1:length(conformer_order)
        if ~options.selected || min(abs(entity.selected-conformer_order(kconf))) == 0
            save_conf = save_conf + 1;
            if options.selected
                curr_pop = entity.populations(conformer_order(kconf))/sum(entity.populations(entity.selected));
            else
                curr_pop = entity.populations(conformer_order(kconf));
            end
            pdbline = sprintf('REMARK 400   MODEL %9i POPULATION %8.4f',save_conf,curr_pop);
            fprintf(fid,'%s\n',pad(pdbline,80));
        end
    end
end

% write ATOM and HETATM records
chains = fieldnames(entity);
select_all = ~options.selected;
indices = zeros(1,5);
save_conf = 0;
for kconf = 1:length(conformer_order)
    atnum = options.datnum;
    % should all conformers be written?
    selected = select_all | all_conformers;
    % or is this specific one selected?
    if min(abs(entity.selected-conformer_order(kconf))) == 0
        selected = true;
        conformer_selected = true;
    else
        conformer_selected = false;
    end
    % skip this conformer, if it is not selected
    if ~selected
        continue
    end
    % write MODEL line, if more than one conformer
    if length(conformer_order) > 1
        if ~conformer_selected || length(entity.selected) > 1 || length(options.order) > 1
            save_conf = save_conf + 1;
            pdbline = sprintf('MODEL %8i',save_conf);
            fprintf(fid,'%s\n',pad(pdbline,80));
        end
    end
    indices(4) = conformer_order(kconf);
    for kc = 1:length(chains)
        chain = chains{kc};
        if isstrprop(chain(1),'upper') % chain fields start with a capital
            selected = entity.(chain).selected | select_all;
%             if conformer_selected
%                 selected = true;
%             end
            indices(1) = entity.(chain).index;
            residues = fieldnames(entity.(chain));
            residue_sorting = sort_residues(residues);
            cid = chain;
            % replace chain identifier, if requested
            if isfield(options,'chainIDs') && ~isempty(options.chainIDs)
                for kx = 1:replace_chains
                    if strcmp(options.chainIDs{kx,1},cid)
                        cid = options.chainIDs{kx,2};
                        break
                    end
                end
            end
            info.cid = cid;
            if selected % expand selection to residues, atoms, and locations, if chain is selected
                for kr = residue_sorting % expand over all residues
                    residue = residues{kr};
                    if strcmp(residue(1),'R') % these are residue fields
                        indices(2) =  entity.(chain).(residue).index;
                        dresnum = 0;
                        if entity.(chain).index <= length(options.dresnum)
                            dresnum = options.dresnum(entity.(chain).index);
                        end
                        info.resnum = str2double(residue(2:end)) + dresnum;
                        % fprintf(1,'(%s).%i\n',chain,info.resnum);
                        info.resname = entity.(chain).(residue).name;
                        info.atomtype = get_atom_type(entity.(chain).(residue).name,biomer);
                        if options.hetframe && strcmp(info.atomtype,'ATOM')
                            continue;
                        end
                        if options.cleaned && strcmp(info.atomtype,'HETATM')
                            continue;
                        end
                        rot_indices = 1:length(entity.(chain).(residue).locations);
                        % rotamers overrule locations
                        if length(entity.(chain).(residue).populations) > 1
                            rot_indices = 1:length(entity.(chain).(residue).populations);
                        end
                        atoms = fieldnames(entity.(chain).(residue));
                        for ka = 1:length(atoms) % expand over all atoms
                            atom = atoms{ka};
                            if isstrprop(atom(1),'upper') && isfield(entity.(chain).(residue).(atom),'index') % these are atom fields
                                indices(3) =  entity.(chain).(residue).(atom).index;
                                atname = atom;
                                if length(atname) > 2 && strcmp(atname(1:2),'Z_')
                                    atname = atname(3:end);
                                end
                                if atname(end) == '_'
                                    atname(end) = '''';
                                end
                                if length(atname) < 2
                                    atname = [atname ' ']; %#ok<AGROW>
                                end
                                if length(atname) < 3
                                    atname = [atname ' ']; %#ok<AGROW>
                                end
                                if length(atname) < 4
                                    atname = [' ' atname]; %#ok<AGROW>
                                end
                                info.atname = atname;
                                info.element = entity.(chain).(residue).(atom).element;
                                if options.bfactor
                                    info.bfactor = entity.(chain).(residue).(atom).bfactor;
                                else
                                    info.bfactor = 0;
                                end
                                if options.charge
                                    info.charge = entity.(chain).(residue).(atom).charge;
                                else
                                    info.charge = [];
                                end
                                loc_str = entity.(chain).(residue).locations;
                                if (loc_str(1) == ' ' || loc_str(1) == '.') && length(loc_str) > 1
                                    loc_str = loc_str(2:end);
                                end
                                for kl = 1:length(rot_indices) % expand over all rotamers/locations
                                    if length(entity.(chain).(residue).(atom).tab_indices) < kl
                                        continue
                                    end
                                    indices(5) = rot_indices(kl); % set current rotamer index
                                    atom_index = entity.(chain).(residue).(atom).tab_indices(kl);
                                    if length(entity.populations) > 1
                                        if length(entity.(chain).(residue).(atom).tab_indices) < length(entity.populations)
                                            continue
                                        end
                                        atom_index = entity.(chain).(residue).(atom).tab_indices(indices(4));
                                    end
                                    info.loc = loc_str(kl);
                                    atnum = wr_pdb_line(fid,entity,atom_index,info,atnum);
                                end
                            end
                        end
                    end
                end
            else % chain was not selected
                for kr = residue_sorting % test all residues
                    residue = residues{kr};
                    if strcmp(residue(1),'R') % these are residue fields
                        indices(2) =  entity.(chain).(residue).index;
                        info.resnum = str2double(residue(2:end));
                        info.resname = entity.(chain).(residue).name;
                        info.atomtype = get_atom_type(entity.(chain).(residue).name,biomer);
                        if options.hetframe && strcmp(info.atomtype,'ATOM')
                            continue;
                        end
                        if options.cleaned && strcmp(info.atomtype,'HETATM')
                            continue;
                        end
                        rot_indices = 1:length(entity.(chain).(residue).locations);
                        % rotamers overrule locations
                        if length(entity.(chain).(residue).populations) > 1
                            rot_indices = 1:length(entity.(chain).(residue).populations);
                        end
                        atoms = fieldnames(entity.(chain).(residue));
                        if entity.(chain).(residue).selected % residue was selected
                            for ka = 1:length(atoms) % expand selection over all atoms
                                atom = atoms{ka};
                                if isstrprop(atom(1),'upper') % these are atom fields
                                    indices(3) =  entity.(chain).(residue).(atom).index;
                                    atname = atom;
                                    if atname(end) == '_'
                                        atname(end) = '''';
                                    end
                                    if length(atname) < 2
                                        atname = [atname ' ']; %#ok<AGROW>
                                    end
                                    if length(atname) < 3
                                        atname = [atname ' ']; %#ok<AGROW>
                                    end
                                    if length(atname) < 4
                                        atname = [' ' atname]; %#ok<AGROW>
                                    end
                                    info.atname = atname;
                                    info.element = entity.(chain).(residue).(atom).element;
                                    if options.bfactor
                                        info.bfactor = entity.(chain).(residue).(atom).bfactor;
                                    else
                                        info.bfactor = 0;
                                    end
                                    if options.charge
                                        info.charge = entity.(chain).(residue).(atom).charge;
                                    else
                                        info.charge = [];
                                    end
                                    loc_str = entity.(chain).(residue).locations;
                                    if loc_str(1) == ' ' && length(loc_str) > 1
                                        loc_str = loc_str(2:end);
                                    end
                                    for kl = 1:length(rot_indices) % expand over all rotamers/locations
                                        indices(:,5) = rot_indices(kl); % set current rotamer index
                                        atom_index = entity.(chain).(residue).(atom).tab_indices(kl);
                                        if length(entity.populations) > 1
                                            if length(entity.(chain).(residue).(atom).tab_indices) < length(entity.populations)
                                                continue
                                            end
                                            atom_index = entity.(chain).(residue).(atom).tab_indices(indices(4));
                                        end
                                        info.loc = loc_str(kl);
                                        atnum = wr_pdb_line(fid,entity,atom_index,info,atnum);
                                    end
                                end
                            end
                        else % no selection on residue level
                            for ka = 1:length(atoms) % test all atoms for selection
                                atom = atoms{ka};
                                if isstrprop(atom(1),'upper') % these are atom fields
                                    indices(3) =  entity.(chain).(residue).(atom).index;
                                    if entity.(chain).(residue).(atom).selected
                                        atname = atom;
                                        if atname(end) == '_'
                                            atname(end) = '''';
                                        end
                                        if length(atname) < 2
                                            atname = [atname ' ']; %#ok<AGROW>
                                        end
                                        if length(atname) < 3
                                            atname = [atname ' ']; %#ok<AGROW>
                                        end
                                        if length(atname) < 4
                                            atname = [' ' atname]; %#ok<AGROW>
                                        end
                                        info.atname = atname;
                                        info.element = entity.(chain).(residue).(atom).element;
                                        if options.bfactor
                                            info.bfactor = entity.(chain).(residue).(atom).bfactor;
                                        else
                                            info.bfactor = 0;
                                        end
                                        if options.charge
                                            info.charge = entity.(chain).(residue).(atom).charge;
                                        else
                                            info.charge = [];
                                        end
                                        loc_str = entity.(chain).(residue).locations;
                                        if loc_str(1) == ' ' && length(loc_str) > 1
                                            loc_str = loc_str(2:end);
                                        end
                                        for kl = 1:length(rot_indices) % expand over all rotamers/locations
                                            indices(5) = rot_indices(kl); % set current rotamer index
                                            atom_index = entity.(chain).(residue).(atom).tab_indices(kl);
                                            if length(entity.populations) > 1
                                                if length(entity.(chain).(residue).(atom).tab_indices) < length(entity.populations)
                                                    continue
                                                end
                                                atom_index = entity.(chain).(residue).(atom).tab_indices(indices(4));
                                            end
                                            info.loc = loc_str(kl);
                                            atnum = wr_pdb_line(fid,entity,atom_index,info,atnum);
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    % write ENDMDL line, if more than one conformer
    if length(conformer_order) > 1
        if ~conformer_selected || length(entity.selected) > 1 || length(options.order) > 1
            fprintf(fid,'%s\n',pad('ENDMDL',80));
        end
    end
end

if nargout < 2
    fprintf(fid,'%s\n',pad('END',80));
    fclose(fid);
end

function atomtype = get_atom_type(name,biomer)
% Determines atom type (ATOM or HETATM) from a residue name using
% definitions of peptide and nucleic acid monomers in biomer

atomtype = 'HETATM';
if contains(upper(biomer.monomers.aa_tags),name)
    atomtype = 'ATOM';
elseif contains(upper(biomer.monomers.nt_tags),name)
    atomtype = 'ATOM';
elseif contains(upper(biomer.monomers.nt_tags_CYANA),name)    
    atomtype = 'ATOM';
end

function atnum = wr_pdb_line(fid,entity,atom_index,info,atnum)
% writes one atom line to the file with identifier fid
% the (1,5) int16 indices refer to entity.index_array in MMMx:atomic format
% atnum is the current atom number in the file
% info has the fields
% .atomtype   either 'ATOM' or 'HETATM'
% .atname     atom name
% .resname    residue name
% .cid        chain identifier
% .resnum     residue number
% .bfactor    temperature factor, B factor
% .charge     charge

[m,~] = size(entity.xyz);

if ~isempty(atom_index) && atom_index > 0 && atom_index <= m
    atnum = atnum + 1;
    xyz = entity.xyz(atom_index,:);
    occ = double(entity.occupancies(atom_index))/100;
    if occ == 0
        occ = 1;
    end
    if isempty(info.charge)
        chg = ' ';
    elseif info.charge > 0
        chg = sprintf('+%i',info.charge);
    elseif info.charge < 0
        chg = sprintf('-%i',abs(info.charge));
    else
        chg = '0';
    end
    if atnum < 100000
        atnum_str = sprintf('%5i',atnum);
    elseif atnum < 1048576
        atnum_str = dec2hex(atnum);
    else
        atnum_str = '*****';
    end
    pdbline = sprintf('%s%s%5s%s%3s%2s%4i%12.3f%8.3f%8.3f%6.2f%6.2f%12s%2s',...
        pad(info.atomtype,6),atnum_str,info.atname,info.loc,info.resname,info.cid,...
        info.resnum,xyz(1),xyz(2),xyz(3),occ,info.bfactor,info.element,chg);
    fprintf(fid,'%s\n',pdbline);
end


function residue_sorting = sort_residues(residues)
% sorts residue field names in acsending order
% non-residue fields are removed in the sorting vector

sorter = zeros(1,length(residues));
true_residues = 0;
for kr = 1:length(residues)
    residue = residues{kr};
    if strcmp(residue(1),'R') % these are residue fields
        sorter(kr) = str2double(residue(2:end)); % residue number
        true_residues = true_residues + 1;
    else
        sorter(kr) = 1e12;
    end
end
[~,residue_sorting] = sort(sorter);
residue_sorting = residue_sorting(1:true_residues);

% function flag = is_backbone(atom_tag)
% 
% peptide_backbone = {'N','CA','C','O'};
% 
% flag = any(strcmp(peptide_backbone,atom_tag));