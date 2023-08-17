function [atom_indices,complete] = get_selection(entity,paradigm)
%
% GET_SELECTION Provide indices for selection
%
%   [atom_indices,complete] = GET_SELECTION(entity)
%   Extracts the current selection in an entity into index arrays 
%
%   [atom_indices,complete] = GET_SELECTION(entity,paradigm)
%   Ignores selection of conformers, rotamers or location if flag paradigm
%   is true, returns only the first loaction in the first rotamer in that
%   case
%
% INPUT
% entity    entity in an MMMx format
% paradigm  flag, if true only the first conformer and rotamer/location is
%           returned, regardless of conformation selections, defaults to
%           false
%
% OUTPUT
% atom_indices  atom index vector, column vector
% complete      complete index matrix at atom level 
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke


max_selected = 100000; % estimate for maximum number of selected objects, 
                       % too small number can make the function slow, a too 
                       % large number makes it memory intensive 

% initialize missing/empty input
if ~exist('paradigm','var') || isempty(paradigm)
    paradigm = false;
end

selection = zeros(max_selected,5,'uint16');

chains = fieldnames(entity);
if paradigm % return only first conformer, regardless of selection
    conformers = 1;
else
    conformers = entity.selected.'; % column vector of selected conformers
end
num_conf = length(conformers);
old_selection = zeros(num_conf,5);
old_selection(1:num_conf,4) = conformers;
old_sel = num_conf;
sel_poi = 0; % pointer into selection array
for kc = 1:length(chains)
    chain = chains{kc};
    selected = false;
    if isstrprop(chain(1),'upper') % chain fields start with a capital
        if entity.(chain).selected
            selected = true;
        end
        old_selection(:,1) = entity.(chain).index;
        residues = fieldnames(entity.(chain));
        if selected % expand selection to residues, atoms, and locations, if chain is selected
            for kr = 1:length(residues) % expand over all residues
                residue = residues{kr};
                if strcmp(residue(1),'R') % these are residue fields
                    res_index =  entity.(chain).(residue).index;
                    rot_indices = 1:length(entity.(chain).(residue).locations);
                    % rotamers overrule locations
                    if length(entity.(chain).(residue).populations) > 1
                        rot_indices = 1:length(entity.(chain).(residue).populations);
                    end
                    if paradigm % return only first rotamer/location, regardless of selection
                        rot_indices = 1;
                    end
                    old_selection(:,2) = res_index; % set current residue index
                    atoms = fieldnames(entity.(chain).(residue));
                    for ka = 1:length(atoms) % expand over all atoms
                        atom = atoms{ka};
                        if isstrprop(atom(1),'upper') % these are atom fields
                            at_index =  entity.(chain).(residue).(atom).index;
                            old_selection(:,3) = at_index; % set current atom index
                            for kl = 1:length(rot_indices) % expand over all rotamers/locations
                                old_selection(:,5) = rot_indices(kl); % set current rotamer index
                                selection(sel_poi+1:sel_poi+old_sel,:) = old_selection; % extend selection
                                sel_poi = sel_poi + old_sel;
                            end
                        end
                    end
                end
            end
        else % chain was not selected
            for kr = 1:length(residues) % test all residues
                residue = residues{kr};
                if strcmp(residue(1),'R') % these are residue fields
                    res_index =  entity.(chain).(residue).index;
                    rot_indices = 1:length(entity.(chain).(residue).locations);
                    % rotamers overrule locations
                    if length(entity.(chain).(residue).populations) > 1
                        rot_indices = 1:length(entity.(chain).(residue).populations);
                    end
                    old_selection(:,2) = res_index; % set current residue index
                    atoms = fieldnames(entity.(chain).(residue));                    
                    if entity.(chain).(residue).selected % residue was selected
                        for ka = 1:length(atoms) % expand selection over all atoms
                            atom = atoms{ka};
                            if isstrprop(atom(1),'upper') % these are atom fields
                                at_index =  entity.(chain).(residue).(atom).index;
                                old_selection(:,3) = at_index; % set current atom index
                                for kl = 1:length(rot_indices) % expand over all rotamers/locations
                                    old_selection(:,5) = rot_indices(kl); % set current rotamer index
                                    selection(sel_poi+1:sel_poi+old_sel,:) = old_selection; % extend selection
                                    sel_poi = sel_poi + old_sel;
                                end
                            end
                        end
                    else % no selection on residue level
                        for ka = 1:length(atoms) % test all atoms for selection
                            atom = atoms{ka};
                            if isstrprop(atom(1),'upper') && isfield(entity.(chain).(residue).(atom),'index') % these are atom fields
                                at_index =  entity.(chain).(residue).(atom).index;
                                old_selection(:,3) = at_index; % set current atom index
                                if entity.(chain).(residue).(atom).selected
                                    for kl = 1:length(rot_indices) % expand over all rotamers/locations
                                        old_selection(:,5) = rot_indices(kl); % set current rotamer index
                                        selection(sel_poi+1:sel_poi+old_sel,:) = old_selection; % extend selection
                                        sel_poi = sel_poi + old_sel;
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
complete = selection(1:sel_poi,:);

[m,~] = size(entity.index_array);
all_atom_indices = 1:m;
[~,atom_indices] = ismember(entity.index_array,complete,'rows');
atom_indices = all_atom_indices(atom_indices~=0);

atom_indices = sort(atom_indices);

complete = entity.index_array(atom_indices,:);

if entity.water_selected
    atom_indices = [atom_indices;entity.water];
end


