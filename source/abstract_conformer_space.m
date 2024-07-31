function [coor,coor_3D,measures] = abstract_conformer_space(entities,addresses)

[D,pop,all_Rg,assignment] = drms_matrix(entities,addresses);
Rg = sqrt(sum(pop.*all_Rg.^2)); 

measures.assignment = assignment;
measures.all_Rg = all_Rg;
measures.populations = pop;

coor = cmdscale(D);
Rg_acs = gyration_radius(coor);

[C,p] = size(coor);

D_check = squareform(pdist(coor));
dev = 2*sqrt(2*sum(sum(triu(D_check-D).^2))/(C*(C-1)));

measures.Rg = Rg;
measures.Rg_acs = Rg_acs;
measures.disorder = Rg_acs/Rg;
measures.dimension = p;
measures.error_nD = dev;

[coor_3D,rmsd,all_rmsd] = refined_3D_embedding(D,all_Rg);

measures.convergence = all_rmsd;
measures.error_3D = rmsd;

function [pair_drms,pop,Rg,assignment] = drms_matrix(entities,addresses)
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
% entities      cell of entities
% addresses     cell of addrss strings selected chains and residue ranges, 
%               optional
%
% OUTPUT
% pair_drms     (C,C) double matrix of backbone pair rmsds (Angstrom)
% pop           (1,C) population vector for C conformers
% Rg            vector (1,c) of radii of gyration
% assignment    (1,C) vector of assignment of conformers to entities
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2024: Gunnar Jeschke

% default input

if ~exist('addresses','var')
    addresses = cell(1,length(entities));
    for ent = 1:length(entities)
        addresses{ent} = '';
    end
end

conformers = zeros(1,length(entities));
for ent = 1:length(entities)
    entity = entities{ent};
    conformers(ent) = length(entity.populations);
end

C = sum(conformers);

assignment = zeros(1,C);
pop = zeros(1,C);
Rg = zeros(1,C);
coor = cell(1,C);
bas = 0;
for ent = 1:length(entities)
    assignment(bas+1:bas+conformers(ent)) = ent;
    entity = entities{ent};
    pop(bas+1:bas+conformers(ent)) = entity.populations;
    [chain,range] = split_chain_range(addresses{ent});
    backbones = get_backbones_ensemble(entity,chain,range);
    chains = fieldnames(backbones);    
    % determine size of complete backbone atoms matrix
    n_atoms = 0;
    for c = 1:length(chains)
        xyz =  backbones.(chains{c}).bb{1};
        [n,~] = size(xyz);
        n_atoms = n_atoms + n;
    end
    for conf = 1:conformers(ent)
        coor0 = zeros(n_atoms,3);
        atom_pointer = 0;
        for c = 1:length(chains)
            xyz =  backbones.(chains{c}).bb{conf};
            [n,~] = size(xyz);
            coor0(atom_pointer+1:atom_pointer+n,:) = xyz;
            atom_pointer = atom_pointer + n;
        end
        coor{bas+conf} = coor0;
        Rg(bas+conf) = gyration_radius(coor0);
    end
    bas = bas + conformers(ent);
    if ent > 1 && n_atoms ~= n_atoms0
        error('abstract_conformer_space : Inconsistent atom number in entity %i',ent);
    end
    n_atoms0 = n_atoms;
end

pair_drms = zeros(C);
N2 = n_atoms*(n_atoms-1)/2;
for c1 = 1:C-1
    for c2 = c1+1:C
        dmat1 = squareform(pdist(coor{c1})); % distance matrix for first conformer
        dmat2 = squareform(pdist(coor{c2})); % distance matrix for second conformer
        drms = sum(sum((dmat1-dmat2).^2)); % mean-square deviation for this conformer pair
        pair_drms(c1,c2) = drms;
        pair_drms(c2,c1) = drms;
    end
end
pair_drms = sqrt(pair_drms/N2); % normalized root mean square deviation
pop = pop/sum(pop); % renormalize populations

function [chain,range] = split_chain_range(address)

chain = '';
range = [];
if isempty(address)
    return
end

poia = strfind(address,'(');
poie = strfind(address,')');
if isempty(poie)
    chain = '';
else
    chain = address(poia+1:poie-1);
    address = address(poie+1:end);
end
if strcmp(chain,'*')
    range = [1,1e6];
    return
end
if isempty(address)
    range = [1,1e6];
    return
end

residues = split(address,'-');
if ~isempty(residues{1})
    range(1) = str2double(residues{1});
end
if ~isempty(residues{2})
    range(2) = str2double(residues{2});
end
if length(range) == 1
    range(2) = range(1);
end