function [pair_drms,pop,Rg,exceptions] = pair_drms_matrix(entity,chain,range,entity2,chain2,range2)
%
% PAIR_DRMS_MATRIX Computes the matrix of conformer pair backbone distance
% root-mean square deviations
%
%   [pair_rmsd,exceptions] = PAIR_DRMS_MATRIX(entity,chain,range)
%   For an entity with C conformes, the C x C matrix of backbone pair
%   distance root-,mean square deviations is computed
%   computation can be confined to a chain or to a residue range within the
%   chain 
%
%   this metric is independent of conformer orientation, see:
%   J. Köfinger, B. Rozycki, G. Hummer, DOI: 10.1007/978-1-4939-9608-7_14
%
%   the MMMx implementation slightly differs by considering all backbone
%   atoms, not one bead per residue
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
% pair_drms     (C,C) double matrix of backbone pair rmsds (Angstrom)
% pop           (1,C) population vector for C conformers
% Rg            vecor (1,c) of radii of gyration
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
%
% exceptions occur if chain length or residue identity is inconsistent
% between conformers in any chain
% however, if all biopolymers in the first PDB file are encountered in
% later files and are consistent, additional biopolymer chains in later
% files do not raise an exception

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
pair_drms = [];

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
    pop = [pop;pop2];
else
    compare = false;
end

Rg = zeros(C,1);

coor = cell(1,C);
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
    coor{conf} = coor0;
    Rg(conf) = gyration_radius(xyz);
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
        coor{C1+conf} = coor0;
        Rg(C1+conf) = gyration_radius(xyz);
    end
end

pair_drms = zeros(C);
N2 = n_atoms*(n_atoms-1)/2;
for k1 = 1:C-1
    for k2 = k1+1:C
        dmat1 = coor2dmat(coor{k1}); % distance matrix for first conformer
        dmat2 = coor2dmat(coor{k2}); % distance matrix for second conformer
        drms = sum(sum((dmat1-dmat2).^2)); % mean-square deviation for this conformer pair
        pair_drms(k1,k2) = drms;
        pair_drms(k2,k1) = drms;
    end
end
pair_drms = sqrt(pair_drms/N2); % normalized root mean square deviation