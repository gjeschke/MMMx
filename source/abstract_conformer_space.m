function [coor,coor_3D,measures,D] = abstract_conformer_space(entities,addresses,options)
% ABSTRACT_CONFORMER_SPACE  Determines the abstract conformer space and its 
%                           3D approximation for a set of ensemble structures
%
%    [coor,coor_3D,measures,D] = abstract_conformer_space(entities,addresses,options)
% 
% INPUT
% entities      cell(1,E) of E ensemble entities
% addresses     cell{1,E) cell of address strings that specify the same
%               polymer chain section in the individual entities, optional,
%               if missing, the whole entity is processed
% options       .convergence    Boolean flag that requests convergence
%                               analysis, defaults to true
%               .sorting        Boolean flag that requests conformer
%                               sorting, defaults to true, can take very
%                               long for large ensembles (1000 or more
%                               conformers)
%               .fname          name of a CSV file for output, defaults to
%                               acs.csv, an empty .fname suppresses file
%                               output, the file is not written if 3D
%                               embedding failed
%               .visualization  name of MMM visualization file, if missing
%                               or empty, no visualization script is
%                               written
%               .graphics       string, 'ribbon' or 'snake' visualization,
%                               defaults to 'snake' (coil with weight
%                               encoded by radius)
%               .subensembles   integer, requests clustering to that number
%                               of subsensembles, defaults to empty (no
%                               clustering)
%
% OUTPUT
% coor          double (C,n) coordinates for C conformers in n-dimensional
%               abstract conformer space (ACS), as determined by classical
%               multidimensional scaling, C is the total number of
%               conformers in all entities
% coor_3D       double(C,3) coordinates for C conformers in an approximate
%               3D ACS; empty if 3D embedding failed
% measures      diagnostic measures on the embedding and 3D approximation,
%               struct with fields
%               .assignment     assignment of conformers to the entities
%               .populations    populations of conformers
%               .Rg             radius of gyration for combined ensemble
%               .all_Rg         radii of gyration for individual conformers
%               .Rg_acs         radius of gyration in n-D ACS
%               .disorder       disorder measure Rg_acs/Rg
%               .dimension      dimension of ACS
%               .extension      extension [Å] in 3D ACS, 3rd root of the
%                               volume of the convec hull, empty if 3D
%                               embedding failed
%               .error_nD       rmsd [Å] of the dRMSD matrix for nD embedding
%               .errors_3D      rmsd errors [Å] of the dRMSD of the
%                               conformers to all other conformers, empty
%                               if 3D embedding failed
%               .mode_3D        origin of 3D coordinates
%               .convergence    convergence of ensemble measures with
%                               ensemble size
%               .sorting        indices that sort the conformers by dRMSD
%                               difference, empty if sorting was not
%                               requested
%               .resolution     resolution in ACS [Å], root mean square of
%                               neareast-neighbor distance
%               .subensembles   assignment of the conformers to
%                               subensembles, empty if no clustering was
%                               requested
% D             dRMSD matrix between the conformers
%
% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% (c) G. Jeschke, 2023-2024

if ~exist('addresses','var')
    addresses = {};
end

if ~exist('options','var') || ~isfield(options,'convergence') || isempty(options.convergence)
    options.convergence = true;
end

if ~isfield(options,'fname')
    options.fname = 'acs.csv';
end

if ~isempty(options.fname) && ~contains(options.fname,'.')
    options.fname = strcat(options.fname,'.csv');
end

if ~isfield(options,'visualization')
    options.visualization = '';
end

if ~isfield(options,'sorting')
    options.sorting = true;
end

if ~isfield(options,'subensembles')
    options.subensembles = [];
end

if ~isfield(options,'graphics') || isempty(options.graphics)
    options.graphics = 'snake';
end

description{1,1} = 'Conformer';
description{1,2} = 'Weight';
description{1,3} = 'x';
description{1,4} = 'y';
description{1,5} = 'z';
description{1,6} = 'Uncertainty';
description{1,7} = 'State';
description{1,8} = 'Cluster';
description{1,9} = 'Cluster type';


[D,pop,all_Rg,assignment] = drms_matrix(entities,addresses);
Rg = sqrt(sum(pop.*all_Rg.^2)); 

[C,~] = size(D);
conformer_order = 1:C;
if options.sorting
    [D,indices] = acs_sorting(D);
    conformer_order = indices;
    pop = pop(indices);
    plot_sorted_distance_matrix(D);
else
    indices = [];
end

measures.assignment = assignment;
measures.sorting = indices;
measures.all_Rg = all_Rg;
measures.populations = pop';

coor = cmdscale(D);
[C,n] = size(coor); % determine number of conformers and dimensionality


colors = turbo(C);

[~,score,~,~,explained] = pca(coor);

if ~isempty(options.subensembles) && options.subensembles > 1
    % compute linkage
    Z = linkage(D,'complete');
    % cluster
    indices = cluster(Z,'maxclust',options.subensembles);
    measures.subensembles = indices;
    colors0 = turbo(options.subensembles+2);
    colors0 = 0.75*colors0(2:end-1,:);
    for c = 1:C
        colors(c,:) = colors0(indices(c),:);
    end
else
    measures.subensembles = [];
    if isempty(indices)
        indices = ones(1,C);
        colors = repmat([0,0.5,0.5],C,1);
    else
        colors = turbo(C);
    end
end

figure; hold on;
for c = 1:C
    plot(score(c,1),score(c,2),'.','MarkerSize',14,'Color',colors(c,:));
end
axis equal
xlabel('x [Å]');
ylabel('y [Å]');
title(sprintf('PCA 2 analysis (%4.1f%% of variance)',sum(explained(1:2))));
set(gca,'FontSize',12);

figure; hold on;
for c = 1:C
    plot3(score(c,1),score(c,2),score(c,3),'.','MarkerSize',14,'Color',colors(c,:));
end
axis equal
xlabel('x [Å]');
ylabel('y [Å]');
zlabel('z [Å]');
title(sprintf('PCA 3 analysis (%4.1f%% of variance)',sum(explained(1:3))));
set(gca,'FontSize',12);
view(30,30);

if ~isempty(options.visualization)
    vname = options.visualization;
    if ~contains(vname,'.')
        vname = strcat(vname,'.mmm');
    end
    mk_visualization_script(vname,entities,assignment,indices,options.graphics);
end

D0 = D + 1e20*eye(C);
nearest_neighbors = min(D0);
measures.resolution = sqrt(sum(pop.*nearest_neighbors.^2));

Rg_acs = gyration_radius(coor);

if options.convergence
    measures.convergence = get_convergence(D,pop);
    measures.convergence.Rg_real = nan(1,C);
    for c = 5:C
        measures.convergence.Rg_real(c) = sqrt(sum(pop(1:c).*all_Rg(1:c).^2)/sum(pop(1:c))); 
    end
end

D_check = squareform(pdist(coor));
dev = sqrt(2*sum(sum(triu(D_check-D).^2))/(C*(C-1)));

measures.Rg = Rg;
measures.Rg_acs = Rg_acs;
measures.disorder = Rg_acs/Rg;
measures.dimension = n;
measures.error_nD = dev;
[coor_3D,~,errcode] = refined_3D_embedding(D,all_Rg,pop'); 
switch errcode
    case 0
        measures.mode_3D = 'non-classical MDS';
    case 1
        measures.mode_3D = 'iterative refinement';
    case 2
        measures.mode_3D = 'classical MDS';
end


if ~isempty(options.fname) && ~isempty(coor_3D)
    if ~isempty(measures.subensembles)
        output_data = zeros(C,9);
    else
        output_data = zeros(C,6);
    end
    output_data(:,1) = conformer_order;
    output_data(:,2) = pop;
    output_data(:,3:5) = coor_3D;
    D_check = squareform(pdist(coor_3D));
    graces = tinv(0.975,C-1)*sqrt(sum((D_check-D).^2/(C-1)));
    output_data(:,6) = graces;
    if ~isempty(measures.subensembles)
        output_data(:,7) = ones(C,1);
        output_data(:,8) = measures.subensembles';
        put_csv(options.fname,output_data,description);
    else
        put_csv(options.fname,output_data,description(:,1:6));
    end
end

if ~isempty(coor_3D)
    [~, vol] = convhulln(coor_3D); % compute the n-dimensional convex hull
    measures.extension = vol^(1/3); % extension in nD ACS 
    
    D_approx = squareform(pdist(coor_3D));
    measures.errors_3D = sqrt(sum(triu(D_approx-D).^2/(C-1)));
    measures.error_3D = sqrt(sum(sum(triu(D_approx-D).^2))/(C*(C-1)));
else
    measures.extension = [];
    measures.errors_3D = [];
    measures.error_3D = [];
end


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

function convergence = get_convergence(D,weights)


[C,~] = size(D);

% randomize sequence of conformers
sorting = rand(1,C);
[~,sorting] = sort(sorting);
D = D(sorting,sorting);
weights = weights(sorting);

convergence.maxmin = nan(1,C);
convergence.rmean = nan(1,C);
convergence.stdr = nan(1,C);
convergence.Rg_acs = nan(1,C);
convergence.dimension = nan(1,C);


if C < 5
    return
end

for c = 5:C
    coor = cmdscale(D,c);
    D_check = squareform(pdist(coor));
    convergence.dimension(c) = sqrt(sum(sum(triu(D_check-D).^2))/(C*(C-1)));
    D0 = D(1:c,1:c);
    coor = cmdscale(D0);
    convergence.Rg_acs(c) = gyration_radius(coor,weights(1:c));
    convergence.rmean(c) = sum(sum(triu(D0)))/(c*(c-1));
    D1 = D0 - convergence.rmean(c)*ones(c);
    convergence.stdr(c) = sqrt(sum(sum(triu(D1).^2))/(c*(c-1)));
    D0 = D0 + 1e20*eye(c);
    nearest_neighbors = min(D0);
    pop = weights(1:c)/sum(weights(1:c));
    convergence.resolution(c) = sqrt(sum(pop.*nearest_neighbors.^2));
    convergence.maxmin(c) = max(min(D0));
end

function mk_visualization_script(fname,entities,assignment,indices,graphics)

if isempty(indices)
    indices = 1:length(assignment);
end
if max(indices) < length(assignment)
    clustered = true;
else
    clustered = false;
end
fid = fopen(fname,'wt')';
fprintf(fid,'%% MMMx sorted conformers visualization script\n');
fprintf(fid,'new !\n');
poi = 0;
nc = zeros(1,10000);
for ke = 1:length(entities)
    entity = entities{ke};
    entity.name = pad(sprintf('E%i',ke),4,'right','X');
    entities{ke} = entity;
    pname = sprintf('acs_entity_%i%.pdb',ke);
    entity.selected = 0;
    put_pdb(entity,pname);
    fprintf(fid,'pdbload %s\n',pname);
    nc(poi+1:poi+length(entity.populations)) = 1:length(entity.populations);
    poi = poi + length(entity.populations);
end
nc = nc(1:poi);
cmap = turbo(max(indices)+2);
for kc = 1:length(indices)
    if clustered
        c = kc;
        n = nc(kc);
    else
        c = indices(kc);
        n = nc(indices(kc));
    end
    ke = assignment(c);
    entity = entities{ke};
    pop = entity.populations(n);
    maxpop = max(entity.populations);
    col = cmap(indices(kc)+1,:);
    switch graphics
        case 'snake'
            fprintf(fid,'show [%s]{%i} coil %6.3f\n',entity.name,n,sqrt(pop/maxpop));
            fprintf(fid,'color [%s]{%i} %6.3f %6.3f %6.3f\n',entity.name,n,col);
        case 'ribbon'
            fprintf(fid,'show [%s]{%i} ribbon\n',entity.name,n);
            fprintf(fid,'color [%s]{%i} %6.3f %6.3f %6.3f\n',entity.name,n,col);
            fprintf(fid,'transparency [%s]{%i} %6.3f\n',entity.name,n,sqrt(pop/maxpop));
    end
end
fclose(fid);

function plot_sorted_distance_matrix(D)

figure;
[n1,n2] = size(D);
plot(1,1,'k.');
plot(n1,n2,'k.');

image(D,'CDataMapping','scaled');

curr_axis = gca;
curr_axis.YDir = 'normal';
c = colorbar;
c.Label.String = 'dRMSD [Å]';
axis tight
xlabel('Conformer');
ylabel('Conformer');
axis equal

