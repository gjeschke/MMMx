function pdm = pair_distance_matrix(entity,chain,range)
%
% PAIR_DISTANCE_MATRIX Computes  avector of all backbone atom pair 
% distances (CA for amino acids or C4' for nucleotides) for each conformer
% and stores it in a matrix [C,p], where C is the number of conformers
% root-mean square deviations
%
%   pdm = PAIR_DISTANCE_MATRIX(entity,chain,range)
%   For an entity with C conformes, the C x p matrix of backbone pair
%   distances, with p = R*(R-1)/2, where R is the number of residues
%
% INPUT
% entity    either an entity or a file name specifying the ensemble with
%           possible extensions
%           .ens   ensemble definition file with populations
%           .pdb   PDB file, filename can contain wildcards, all files 
%                  matching it are processed, 
%           none   extension .pdb is appended
% chain     selected chain, optional
% range     range of residues, by number integer(1,2), optional
%
% OUTPUT
% pdm       (C,p) double matrix of pair distances

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2022: Gunnar Jeschke

% default input

if ~exist('chain','var')
    chain = '';
end

if ~exist('range','var')
    range = [];
end

% initialize empty output
pdm = [];

[backbones,~,exceptions] = get_backbones_ensemble(entity,chain,range);

if ~isempty(exceptions{1})
    return
end

C = length(entity.populations);

coor = cell(1,C);
chains = fieldnames(backbones);

% determine size of complete backbone atoms matrix
n_atoms = 0;
for c = 1:length(chains)
    xyz =  backbones.(chains{c}).bb{1};
    [n,~] = size(xyz);
    n_atoms = n_atoms + n;
end

for conf = 1:C
    coor0 = zeros(n_atoms,3);
    atom_pointer = 0;
    for c = 1:length(chains)
        xyz =  backbones.(chains{c}).bb{conf};
        [n,~] = size(xyz);
        coor0(atom_pointer+1:atom_pointer+n,:) = xyz;
        atom_pointer = atom_pointer + n;
    end
    coor{conf} = coor0;
end

N2 = n_atoms*(n_atoms-1)/2;
pdm = zeros(C,N2);
for c = 1:C
    dmat = coor2dmat(coor{c}); % distance matrix for first conformer
    D = diag(dmat);
    pdm(c,:) = squareform((dmat-diag(D)).');
end