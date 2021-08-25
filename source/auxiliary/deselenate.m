function [entity,num_Se] = deselenate(entity)
%  [entity,num_Se] = DESELENATE(entity)
%    Replaces selenocysteine and selenomethionine by cysteine and
%    methionine, respectively
%
% INPUT:
%
% entity        input entity
%
% OUTPUT:
%
% entity        entity with replaced sidechains
% num_Se        number of replace seleno sidechains
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

num_Se = 0;
chains = fieldnames(entity);
for ch = 1:length(chains)
    chain = chains{ch};
    if isstrprop(chain(1),'upper') % chain fields start with a capital
        residues = fieldnames(entity.(chain));
        for kr = 1:length(residues) % expand over all residues
            residue = residues{kr};
            if strcmp(residue(1),'R') % these are residue fields
                if strcmpi(entity.(chain).(residue).name,'CSE')
                    entity.(chain).(residue).SG = entity.(chain).(residue).SE;
                    entity.(chain).(residue).name = 'CYS';
                    entity.(chain).(residue) = rmfield(entity.(chain).(residue),'SE');
                    entity = correct_pos_elm(entity,entity.(chain).(residue).SG);
                    num_Se = num_Se + 1;
                end
                if strcmpi(entity.(chain).(residue).name,'MSE')
                    entity.(chain).(residue).SD = entity.(chain).(residue).SE;
                    entity.(chain).(residue).name = 'MET';
                    entity.(chain).(residue) = rmfield(entity.(chain).(residue),'SE');
                    entity = correct_pos_elm(entity,entity.(chain).(residue).SD);
                    num_Se = num_Se + 1;
                end
            end
        end
    end
end

function entity = correct_pos_elm(entity,atom)
% corrects position and element number from Se to S

max_bond = 2; % maximum length of a C-Se bond

for k = 1:length(atom.tab_indices)
    at = atom.tab_indices(k);
    entity.elements(at) = 16;
    dist = sqrt(sum((entity.xyz - entity.xyz(at,:)).^2,2));
    indices = find(dist < max_bond);
    vec = -0.1*(entity.xyz(at,:) - entity.xyz(indices(1),:)); % S-C bond is 10% shorter than Se-C bond
    vec2 = -0.1*(entity.xyz(indices(3),:) - entity.xyz(at,:)); % S-C bond is 10% shorter than Se-C bond
    entity.xyz(at,:) = entity.xyz(at,:) + vec;
    entity.xyz(indices(3),:) = entity.xyz(indices(3),:) + vec + vec2;
end
    
