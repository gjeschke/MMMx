function [entity,eau] = ensemble_aligned_uncertainty(entity)
%
% ENSEMBLE_ALIGNED_UNCERTAINTY   Computes the ensemble aligned uncertainty
%                                matrix for an entity and adds it to the
%                                entity
%
%   [entity,eau] = ensemble_aligned_uncertainty(entity)
%
%   WARNING:    this is currently implemented only for entities consisting
%               of a single protein chain
%
% INPUT
% entity        MMMx:atomic entity; the current version supports only
%               single-chain entities
%
% OUTPUT
% entity       entity structure in MMMx:atomic representation 
%              the input entity is augmented by
%              .eau             EAU matrix
%              .eau_resnums     residue numbers for EAU matrix
%              .eau_chain       chain identifier
% eau          EAU matrix
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2025: Gunnar Jeschke

% analyze the first chain
chains = fieldnames(entity);
for c = 1:length(chains)
    chain = chains{c};
    if isstrprop(chain(1),'upper')
        break
    end
end

residues = fieldnames(entity.(chain));
resmin = 1e6;
resmax = 0;
for kr = 1:length(residues)
    residue = residues{kr};
    if strcmp(residue(1),'R') % these are residue fields
        resnum = str2double(residue(2:end));
        if resnum < resmin
            resmin = resnum;
        end
        if resnum > resmax
            resmax = resnum;
        end
    end
end

resnums = resmin:resmax;
eau = zeros(length(resnums));

for refres = resnums
    selected = sprintf('(%s)%i',chain,refres);
    centity = superimpose_ensemble(entity,selected);
    for distributed = resnums
        residue = sprintf('R%i',distributed);
        indices = centity.(chain).(residue).CA.tab_indices;
        coor = centity.xyz(indices,:);
        coor = coor - repmat(mean(coor),length(indices),1);
        rmsd = sqrt(sum(sum(coor.^2))/length(indices));
        eau(refres,distributed) = 2*rmsd;
    end
end

entity.eau = eau(resnums,resnums);
entity.eau_resnums = resnums;
entity.eau_chain = chain;