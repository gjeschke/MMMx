function [coor,coor_3D,measures,D] = abstract_conformer_space(entities,addresses)
% ABSTRACT_CONFORMER_SPACE  Determines the abstract conformer space and its 
%                           3D approximation for a set of ensemble structures
%
%    [coor,coor_3D,measures] = abstract_conformer_space(entities,addresses)
% 
% INPUT
% entities      cell(1,E) of E ensemble entities
% addresses     cell{1,E) cell of address strings that specify the same
%               polymer chain section in the individual entities, optional,
%               if missing, the whole entity is processed
%
% OUTPUT
% coor          double (C,n) coordinates for C conformers in n-dimensional
%               abstract conformer space (ACS), as determined by classical
%               multidimensional scaling, C is the total number of
%               conformers in all entities
% coor_3D       double(C,3) coordinates for C conformers in an approximate
%               3D ACS
% measures      diagnostic measures on the embedding and 3D approximation,
%               struct with fields
%               .assignment     assignment of conformers to the entities
%               .populations    populations of conformers
%               .Rg             radius of gyration for combined ensemble
%               .all_Rg         radii of gyration for individual conformers
%               .Rg_acs         radius of gyration in n-D ACS
%               .disorder       disorder measure Rg_acs/Rg
%               .dimension      dimension of ACS
%               .extension      extension [Å] in nD ACS, n-th root of the
%                               hypervolume of the convec hull
%               .error_nD       rmsd [Å] of the dRMSD matrix for nD embedding
%               .errors_3D      rmsd errors [Å] of the dRMSD of the
%                               conformers to all other conformers
%               .error_3D       rmsd [Å] of the dRMSD matrix for 3D embedding
%               .convergence    convergence of ensemble measures with
%                               ensemble size
% D             dRMSD matrix between the conformers
%
% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% (c) G. Jeschke, 2023-2024

if ~exist('addresses','var')
    addresses = {};
end

[D,pop,all_Rg,assignment] = drms_matrix(entities,addresses);
Rg = sqrt(sum(pop.*all_Rg.^2)); 

measures.assignment = assignment;
measures.all_Rg = all_Rg;
measures.populations = pop';

coor = cmdscale(D);
[C,n] = size(coor); % determine number of conformers and dimensionality
measures.convergence = get_convergence(D);

Rg_acs = gyration_radius(coor);
figure(1); clf; hold on;
plot(1:C,measures.convergence.rmean,'k.','MarkerSize',12);
plot(1:C,measures.convergence.stdr,'b.','MarkerSize',12);
plot(1:C,measures.convergence.maxmin,'r.','MarkerSize',12);

D_check = squareform(pdist(coor));
dev = sqrt(2*sum(sum(triu(D_check-D).^2))/(C*(C-1)));

measures.Rg = Rg;
measures.Rg_acs = Rg_acs;
measures.disorder = Rg_acs/Rg;
measures.dimension = n;
measures.error_nD = dev;
coor_3D = refined_3D_embedding(D,all_Rg,pop'); 

[~, vol] = convhulln(coor_3D); % compute the n-dimensional convex hull
measures.extension = vol^(1/3); % extension in nD ACS 

D_approx = squareform(pdist(coor_3D));
measures.errors_3D = sqrt(sum(triu(D_approx-D).^2/(C-1)));


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

if ~exist('addresses','var') || isempty(addresses)
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

function convergence = get_convergence(D)

[C,~] = size(D);
convergence.maxmin = nan(1,C);
convergence.rmean = nan(1,C);
convergence.stdr = nan(1,C);

if C < 5
    return
end

for c = 5:C
    D0 = D(1:c,1:c); 
    convergence.rmean(c) = sum(sum(triu(D0)))/(c*(c-1));
    D1 = D0 - convergence.rmean(c)*ones(c);
    convergence.stdr(c) = sqrt(sum(sum(triu(D1).^2))/(c*(c-1)));
    D0 = D0 + 1e20*eye(c);
    convergence.maxmin(c) = max(min(D0));
end


