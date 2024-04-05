function [pair_rmsd,pop,exceptions,Rg] = pair_rmsd_matrix(entity,chain,range,entity2,chain2,range2)
%
% PAIR_RMSD_MATRIX Computes the matrix of conformer pair backbone root mean
%                   square deviations upon optimal superposition
%
%   [pair_rmsd,exceptions] = PAIR_RMSD_MATRIX(entity,chain,range)
%   For an entity with C conformes, the C x C matrix of backbone pair rmsds
%   is computed
%   computation can be confined to a chain or to a residue range within the
%   chain 
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
% The following inputs can be used for comparing two ensembles
%
% entity2   either an entity or a file name specifying the ensemble with
%           possible extensions
%           .ens   ensemble definition file with populations
%           .pdb   PDB file, filename can contain wildcards, all files 
%                  matching it are processed, 
%           none   extension .pdb is appended
% chain2    selected chain, optional
% range2    range of residues, by number integer(1,2), optional
%
% OUTPUT
% pair_rmsd     (C,C) double matrix of backbone pair rmsds (Angstrom)
% pop           (1,C) population vector for C conformers
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
% Rg            (1,C) vector of the radii of gyration for C conformers 
%
% exceptions occur if chain length or residue identity is inconsistent
% between conformers in any chain
% however, if all biopolymers in the first PDB file are encountered in
% later files and are consistent, additional biopolymer chains in later
% files do not raise an exception

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% default input

Rg = [];

if ~exist('chain','var')
    chain = '';
end

if ~exist('range','var')
    range = [];
end

% initialize empty output
pair_rmsd = [];

[backbones,pop,exceptions] = get_backbones_ensemble(entity,chain,range);

if ~isempty(exceptions{1})
    return
end

C1 = length(pop);
C = C1;

if exist('entity2','var') && ~isempty(entity2) % ensemble comparison    
    compare = true;    
    if ~exist('chain2','var')
        chain2 = '';
    end    
    if ~exist('range2','var')
        range2 = [];
    end    
    [backbones2,pop2,exceptions] = get_backbones_ensemble(entity2,chain2,range2);    
    if ~isempty(exceptions{1})
        return
    end
    C2 = length(pop2);
    C = C1 + C2;
else
    compare = false;
end

coor = cell(1,C);
Rg = zeros(C,1);
chains = fieldnames(backbones);

% determine size of complete backbone atoms matrix
n_atoms = 0;
for c = 1:length(chains)
    xyz =  backbones.(chains{c}).bb{1};
    [n,~] = size(xyz);
    n_atoms = n_atoms + n;
end

if compare
    chains2 = fieldnames(backbones2);
    for c = 1:length(chains2)
        xyz =  backbones2.(chains2{c}).bb{1};
        [n,~] = size(xyz);
        n_atoms = n_atoms + n;
    end
end

for conf = 1:C1
    coor0 = zeros(n_atoms,3);
    atom_pointer = 0;
    for c = 1:length(chains)
        xyz =  backbones.(chains{c}).bb{conf};
        [n,~] = size(xyz);
        coor0(atom_pointer+1:atom_pointer+n,:) = xyz;
        atom_pointer = atom_pointer + n;
    end
    Rg(conf) = gyration_radius(coor0);
    coor{conf} = coor0;
end

if compare
    for conf = 1:C2
        coor0 = zeros(n_atoms,3);
        atom_pointer = 0;
        for c = 1:length(chains2)
            xyz =  backbones2.(chains2{c}).bb{conf};
            [n,~] = size(xyz);
            coor0(atom_pointer+1:atom_pointer+n,:) = xyz;
            atom_pointer = atom_pointer + n;
        end
        Rg(C1+conf) = gyration_radius(coor0);
        coor{C1+conf} = coor0;
    end
end

pair_rmsd = zeros(C);
for k1 = 1:C-1
    for k2 = k1+1:C
        rmsd = rmsd_superimpose(coor{k1},coor{k2});
        pair_rmsd(k1,k2) = rmsd;
        pair_rmsd(k2,k1) = rmsd;
    end
end
