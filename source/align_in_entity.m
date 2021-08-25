function res_corr = align_in_entity(entity1,chain1,entity2,chain2)
%
% ALIGN_IN_ENTITY Align two peptide chains in two entities
%
%   residue_correspondence = align_in_entity(entity1,chain1,entity2,chain2)
%   Finds pairs of residues that are aligned between chain1 in entity1 and
%   chain2 in entity2, uses MUSCLE by R. C. Edgar for alignment
%
% INPUT:
% 
% entity1   first entity in MMMx:atomic format
% chain1    chain in first entity
% entity2   second entity in MMMx:atomic format
% chain2    chain in second entity
%
% OUTPUT:
%
% res_corr  double(r,2), r residue number pairs of aligned residues
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% initialize empty output
res_corr = [];

% determine the two sequences
residues1 = fieldnames(entity1.(chain1));
sequence1 = char(zeros(1,length(residues1)));
resnum1 = zeros(1,length(residues1));
spoi = 0;
for kr = 1:length(residues1)
    residue = residues1{kr};
    if strcmp(residue(1),'R') % these are residue fields
        slc = tlc2slc(entity1.(chain1).(residue).name);
        if ~isempty(slc)
            spoi = spoi + 1;
            sequence1(spoi) = slc;
            resnum1(spoi) = str2double(residue(2:end));
        end
    end
end
sequence1 = sequence1(1:spoi);
resnum1 = resnum1(1:spoi);

residues2 = fieldnames(entity2.(chain2));
sequence2 = char(zeros(1,length(residues2)));
resnum2 = zeros(1,length(residues2));
spoi = 0;
for kr = 1:length(residues2)
    residue = residues2{kr};
    if strcmp(residue(1),'R') % these are residue fields
        slc = tlc2slc(entity2.(chain2).(residue).name);
        if ~isempty(slc)
            spoi = spoi + 1;
            sequence2(spoi) = slc;
            resnum2(spoi) = str2double(residue(2:end));
        end
    end
end
sequence2 = sequence2(1:spoi);
resnum2 = resnum2(1:spoi);

seqs{1} = sequence1;
seqs{2} = sequence2;
alignment = sequence_alignment(seqs);
[~,N] = size(alignment);
res_corr = zeros(N,2);
poi = 0;
for k = 1:N
    if ~isnan(alignment(1,k)) && ~isnan(alignment(2,k)) % this is an aligned residue pair
        poi = poi + 1;
        res_corr(poi,1) = resnum1(alignment(1,k));
        res_corr(poi,2) = resnum2(alignment(2,k));
    end
end
res_corr = res_corr(1:poi,:);
