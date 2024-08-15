function backbone = get_backbone(entity,chain,conformer,range,options)
%
% GET_BACKBONE Retrieve backbone coordinates for a chain (one conformer)
%
%   backbone = GET_BACKBONE(entity,chain,conformer,range,options)
%   Provides coordinates and corresponding residue numbers for all CA and
%   C4' atoms in one conformer of a chain, optionally, a residue range can
%   be specified
%
%
% INPUT
% entity    entity in an MMMx format
% chain     chain identifier, defaults to A, it is discouraged to not
%           supply this argument
% conformer conformer (model) number, defaults to 1
% range     double(1,2), optional residue range, defaults to all residues
% options   .full   flag, if true, N and C coordinates for amino acid
%                   residues are also extracted
%           .list   if true, range argument is interpreted as a list of residues
%                   defaults to false
%
% OUTPUT
% backbone      struct with fields, all fields are empty if the chain is
%               not a biopolymer
%               .aa  (1,naa) double, numbers of residues with a CA atom
%               .CA  (naa,3) coordinate array for CA atoms
%               .res string with length naa, amino acid single-letter codes
%               .nt  (1,nnt) double, numbers of residues with a C4' atom
%               .C4p (nnt,3) coordinate array for C4' atoms
%               .nuc string with length ntt, nucleotide single-letter codes
%               if options.full is true
%               .N   (naa,3) coordinate array for N atoms
%               .C   (naa,3) coordinate array for C atoms
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% load monomer definitions
defs = load('monomers.mat');

if ~exist('options','var') || isempty(options) || ~isfield(options,'full')
    options.full = false;
end
if ~isfield(options,'list')
    options.list = false;
end

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

if ~exist('range','var') || isempty(range)
    range = [-1e10,1e10];
end

maxres = 10000; % maximum number of residues
aa = zeros(1,maxres);
CA = zeros(maxres,3);
if options.full
    N = zeros(maxres,3);
    C = zeros(maxres,3);
end
res = blanks(maxres);
aa_poi = 0; % counter for detected CA atoms
nt = zeros(1,maxres);
C4p = zeros(maxres,3);
nuc = blanks(maxres);
nt_poi = 0; % counter for detected C4' atoms

residues = fieldnames(entity.(chain));
for kr = 1:length(residues) % loop over all residues
    residue = residues{kr};
    if strcmp(residue(1),'R') % these are residue fields
        resnum = str2double(residue(2:end)); % residue number
		selected = false;
		if options.list
		    if min(abs(range-resnum)) == 0
			    selected = true;
		    end
		elseif resnum >= range(1) && resnum <= range(2)
		    selected = true;
		end
        if selected % is the residue in the requested range
            tlc = entity.(chain).(residue).name;
            if isfield(entity.(chain).(residue),'CA') % is there a CA atom?
                aa_poi = aa_poi + 1;
                aa(aa_poi) = resnum;
                CA(aa_poi,:) = entity.xyz(entity.(chain).(residue).CA.tab_indices(conformer),:);
                if isfield(defs.residues,tlc)
                    res(aa_poi) = defs.residues.(tlc).slc;
                else
                    res(aa_poi) = '?';
                end
                if options.full
                    N(aa_poi,:) = entity.xyz(entity.(chain).(residue).N.tab_indices(conformer),:);
                    C(aa_poi,:) = entity.xyz(entity.(chain).(residue).C.tab_indices(conformer),:);
                end
            end
            if isfield(entity.(chain).(residue),'C4_') % is there a C4' atom?
                nt_poi = nt_poi + 1;
                nt(nt_poi) = resnum;
                C4p(nt_poi,:) = entity.xyz(entity.(chain).(residue).C4_.tab_indices(conformer),:);
                if isfield(defs.residues,tlc)
                    nuc(nt_poi) = defs.residues.(tlc).slc;
                else
                    nuc(nt_poi) = '?';
                end
            end
        end
    end
end

% assign output
backbone.aa = aa(1:aa_poi);
backbone.CA = CA(1:aa_poi,:);
backbone.res = res(1:aa_poi);
backbone.nt = nt(1:nt_poi);
backbone.C4p = C4p(1:nt_poi,:);
backbone.nuc = nuc(1:nt_poi);

if options.full
    backbone.N = N(1:aa_poi,:);
    backbone.C = C(1:aa_poi,:);
end
