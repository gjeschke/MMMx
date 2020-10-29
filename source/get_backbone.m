function backbone = get_backbone(entity,chain,conformer)
%
% GET_BACKBONE Retrieve backbone coordinates for a chain (one conformer)
%
%   backbone = GET_BACKBONE(entity,chain,conformer)
%   Provides coordinates and corresponding residue numbers for all CA and
%   C4' atoms in one conformer of a chain
%
%
% INPUT
% entity    entity in an MMMx format
% chain     chain identifier, defaults to A, it is discouraged to not
%           supply this argument
% conformer conformer (model) number, defaults to 1
%
% OUTPUT
% backbone      struct with fields, all fields are empty if the chain is
%               not a biopolymer
%               .aa  (1,naa) double, numbers of residues with a CA atom
%               .CA  (naa,3) coordinate array for CA atoms
%               .nt  (1,nnt) double, numbers of residues with a C4' atom
%               .C4p (nnt,3) coordinate array for C4' atoms
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
backbone.aa = [];
backbone.CA = [];
backbone.nt = [];
backbone.C4p = [];

% set default input
if ~exist('chain','var') || isempty(chain)
    chain = 'A';
end

if ~exist('conformer','var') || isempty(conformer)
    conformer = 1;
end

maxres = 10000; % maximum number of residues
aa = zeros(1,maxres);
CA = zeros(maxres,3);
aa_poi = 0; % counter for detected CA atoms
nt = zeros(1,maxres);
C4p = zeros(maxres,3);
nt_poi = 0; % counter for detected C4' atoms

residues = fieldnames(entity.(chain));
for kr = 1:length(residues) % loop over all residues
    residue = residues{kr};
    if strcmp(residue(1),'R') % these are residue fields
        resnum = str2double(residue(2:end)); % residue number
        if isfield(entity.(chain).(residue),'CA') % is there a CA atom?
            aa_poi = aa_poi + 1;
            aa(aa_poi) = resnum;
            CA(aa_poi,:) = entity.xyz(entity.(chain).(residue).CA.tab_indices(conformer),:);
        end
        if isfield(entity.(chain).(residue),'C4_') % is there a C4' atom?
            nt_poi = nt_poi + 1;
            nt(nt_poi) = resnum;
            C4p(nt_poi,:) = entity.xyz(entity.(chain).(residue).C4_.tab_indices(conformer),:);
        end
    end
end

% assign output
backbone.aa = aa(1:aa_poi);
backbone.CA = CA(1:aa_poi,:);
backbone.nt = nt(1:nt_poi);
backbone.C4p = C4p(1:nt_poi,:);
