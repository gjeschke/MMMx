function [atoms,objects] = get_selection(entity)
%
% GET_SELECTION Provide indices for selection (MMMx internal) 
%
%   [atoms,objects] = GET_SELECTION
%   Extracts the current selection in an entity into index arrays 
%
%   [atoms,objects] = GET_SELECTION(entity)
%   Scans an entity for selected objects on all hierachy levels and 
%   builds an atom index vector as well as an object index matrix
%
% INPUT
% entity    entity in an MMMx format
%
% OUTPUT
% atoms     atom index vector, column vector
% objects   object index matrix 
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

profile on

max_selected = 100000; % estimate for maximum number of selected objects, 
                       % too small number can make the function slow, a too 
                       % large number makes it memory intensive 

atoms = [];
objects = zeros(max_selected,5);
selection = zeros(max_selected,5);
complete = [];

chains = fieldnames(entity);
conformers = entity.selected.'; % column vector of selected conformers
num_conf = length(conformers);
sel_poi = 0; % pointer into selection array
for kc = 1:length(chains)
    chain = chains{kc};
    selected = false;
    if isstrprop(chain(1),'upper') % chain fields start with a capital
        selection(1:num_conf,1) =  entity.(chain).index;
        selection(1:num_conf,4) = conformers;
        sel_poi = num_conf;
        if entity.(chain).selected
            selected = true;
        end
        residues = fieldnames(entity.(chain));
        old_selection = selection(1:sel_poi,:);
        old_sel = sel_poi;
        if selected % expand selection, if chain is selected
            for kr = 1:length(residues) % expand over all residues
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
                    if entity.(chain).(residue).selected
                        res_index =  entity.(chain).(residue).index;
                        rot_indices = 1:length(entity.(chain).(residue).locations);
                    end
                end
            end
        end
    end
end

if entity.water_selected
    atoms = [atoms;entity.water];
end

profile viewer
