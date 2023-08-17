function clusters = cluster_transition(entity1,entity2,options)
%
% CLUSTER_TRANSITION Characterize and visualize transition between ensemble
%                       structures via clustering of conformers
%
%   clusters = CLUSTER_TRANSITION(entity1,entity2,options)
%   Clusters all conformers of two ensembles, assigns clusters as belonging
%   to ensemble 1 or ensemble 2 or mixed. Splits the mixed clusters into
%   subclusters. Determines similarities for all cluster and subcluster
%   pairs. Writes a visualization script for MMM if requested.
%
% INPUT
% entity1       first input entity, must be an ensemble with several
%               conformers
% entity2       second input entity, must be an ensemble with several
%               conformers and must have the same primary structure as
%               entity 1
% options       options, struct
%               .nc             number of clusters, defaults to 9
%               .visualize      name of visualization script, if empty, no
%                               visualization script is written
%               .superimposed   range of superimposed residues, needed only
%                               for visualization script, defaults to empty
%                               string
%               .chain1         chain selection for entity 1, defaults to
%                               empty string (all chains)
%               .chain2         chain selection for entity 2, defaults to
%                               empty string (all chains)
%               .fname1         basis name for PDB file of ensemble 1, used
%                               only if a visualization script is written
%               .fname2         basis name for PDB file of ensemble 2, used
%                               only if a visualization script is written
% OUTPUT
% clusters      information on clusters, struct
%               .nc             number of clusters
%               .C1             number of conformers in first ensemble
%               .C2             number of conformers in second ensemble
%               .members        cell(nc,2), members of each cluster
%                               belonging to ensemble 1 and to ensemble 2
%               .type           cluster type, 1 for pure ensemble 1, 2 for
%                               pure ensemble 2, 0 for mixed
%               .similarities   double matrix (nc,nc) of similarity
%                               measures
%               .subsimilar     double vector (nc,1) of similarities
%                               between subclusters, NaN for pure clusters
%                               or if a cluster contains less then three
%                               conformers of one of the two ensembles
%               .scale          double vector (nc,1) of a one-dimensional
%                               similarity scale between 0 and 1, with 0
%                               corresponding to "most" ensemble-1-like
%                               cluster and 1 to "most" ensemble-2-like
%                               cluster, mixed clusters fall between 1/3
%                               and 2/3
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

if ~exist('options','var') || isempty(options) || ~isfield(options,'nc')
    options.nc = 9;
end

if ~isfield(options,'superimposed')
    options.superimposed = '';
end

if ~isfield(options,'chain1')
    options.chain1 = '';
end

if ~isfield(options,'chain2')
    options.chain2 = '';
end

clusters.C1 = length(entity1.populations);
clusters.C2 = length(entity2.populations);
populations = [entity1.populations;entity2.populations];
D = pair_drms_matrix(entity1,options.chain1,[],entity2,options.chain2,[]);
% compute linkage
Z = linkage(D,'average');
valid = false;
while ~valid
    % cluster
    assignment = cluster(Z,'maxclust',options.nc);
    valid = true;
%     for k = 1:options.nc
%         if sum(assignment==k) < 3
%             valid = false;
%         end
%     end
    if ~valid
        options.nc = options.nc - 1;
    end
end
clusters.nc = options.nc;
clusters.members = cell(clusters.nc,2);
clusters.type = zeros(clusters.nc,1);


% assign clusters, cluster types 1, 2 and 0 (mixed)
conformers = 1:length(assignment);
for clust = 1:clusters.nc
    members = conformers(assignment == clust);
    if max(members) <= clusters.C1 % pure ensemble 1 cluster
        clusters.type(clust) = 1;
        clusters.members{clust,1} = members;
        clusters.members{clust,2} = [];
    end
    if min(members) > clusters.C1 % pure ensemble 2 cluster
        clusters.type(clust) = 2;
        clusters.members{clust,1} = [];
        clusters.members{clust,2} = members - clusters.C1;
    end
    if clusters.type(clust) == 0 % mixed cluster
        clusters.members{clust,1} = members(members <= clusters.C1);
        clusters.members{clust,2} = members(members > clusters.C1) - clusters.C1;
    end
end

clusters.similarities = ones(clusters.nc);
for cl1 = 1:clusters.nc-1
    indices1 = [clusters.members{cl1,1} clusters.members{cl1,2} + clusters.C1];
    for cl2 = cl1+1:clusters.nc
        indices2 = [clusters.members{cl2,1} clusters.members{cl2,2} + clusters.C1];
        clusters.similarities(cl1,cl2) = get_similarity(D,populations,indices1,indices2);
        clusters.similarities(cl2,cl1) = clusters.similarities(cl1,cl2);
    end
end

clusters.subsimilar = NaN(clusters.nc,1);
for clust = 1:clusters.nc
    if clusters.type(clust) == 0
        indices1 = clusters.members{clust,1};
        indices2 = clusters.members{clust,2} + clusters.C1;
        if length(indices1) >= 3 && length(indices2) >=3
            clusters.subsimilar(clust) = get_similarity(D,populations,indices1,indices2);
        end
    end
end

scaled = zeros(clusters.nc,1);
for cl1 = 1:clusters.nc
    for cl2 = 1:clusters.nc
        if clusters.type(cl2) == 1
            scaled(cl1) = scaled(cl1) - clusters.similarities(cl1,cl2);
        end
        if clusters.type(cl2) == 2
            scaled(cl1) = scaled(cl1) + clusters.similarities(cl1,cl2);
        end
    end
end

% transform to scale bewteen 0 and 1
scaled = scaled - min(scaled);
scaled = scaled/max(scaled);

rescaled = scaled;
for clust = 1:clusters.nc
    type = clusters.type(clust);
    x0 = scaled(clusters.type == type);
    x(1) = min(x0);
    x(2) = max(x0);
    switch type
        case 0
            if x(1) == x(2)
                rescaled(clust) = 0.5;
            else
                rescaled(clust) = 0.4+0.2*(scaled(clust)-x(1))/(x(2)-x(1));
            end
        case 1
            if x(1) == x(2)
                rescaled(clust) = 0;
            else
                rescaled(clust) = 1/5*(scaled(clust)-x(1))/(x(2)-x(1));
            end
        case 2
            if x(1) == x(2)
                rescaled(clust) = 1;
            else
                rescaled(clust) = 1/5*(4+(scaled(clust)-x(1))/(x(2)-x(1)));
            end
    end
end

colors = cluster_color_map;


figure(1); clf; hold on;

for clust = 1:clusters.nc
    switch clusters.type(clust)
        case 0
            marker = 'o';
        case 1
            marker = 'v';
        case 2
            marker = '^' ;
    end
    color_index = 1 + round(100*rescaled(clust));
    col = colors(color_index,:);
    plot(scaled(clust),rescaled(clust),marker,'MarkerSize',10,'Color',col,'MarkerFaceColor',col);
end

if isfield(options,'visualize') && ~isempty(options.visualize)
    [pname,bname,ext] = fileparts(options.fname1);
    if isempty(ext)
        options.fname1 = fullfile(pname,sprintf('%s.pdb',bname));
    end
    [pname,bname,ext] = fileparts(options.fname2);
    if isempty(ext)
        options.fname2 = fullfile(pname,sprintf('%s.pdb',bname));
    end
    entity1.name = 'DISO';
    entity2.name = 'ORDE';
    put_pdb(entity1,options.fname1);
    put_pdb(entity2,options.fname2);
    [pname,bname,ext] = fileparts(options.visualize);
    graphics_x =  fullfile(pname,sprintf('%s_view_x.png',bname));
    graphics_y =  fullfile(pname,sprintf('%s_view_y.png',bname));
    if isempty(ext)
        options.visualize = fullfile(pname,sprintf('%s.mmm',bname));
    end
    ofid = fopen(options.visualize,'wt');
    fprintf(ofid,'%% MMMx transition visualization script\n');
    fprintf(ofid,'new !\n');
    fprintf(ofid,'pdbload %s\n',options.fname1);
    fprintf(ofid,'pdbload %s\n',options.fname2);
    for clust = 1:clusters.nc
        color_index = 1 + round(100*rescaled(clust));
        col = colors(color_index,:);
        conformers = clusters.members{clust,1};
        for c = 1:length(conformers)
            pop = entity1.populations(conformers(c));
            fprintf(ofid,'show [DISO]{%i} coil %6.3f\n',conformers(c),sqrt(pop));
            fprintf(ofid,'color [DISO]{%i} %6.3f %6.3f %6.3f\n',conformers(c),col);
        end
        conformers = clusters.members{clust,2};
        for c = 1:length(conformers)
            pop = entity2.populations(conformers(c));
            fprintf(ofid,'show [ORDE]{%i} coil %6.3f\n',conformers(c),sqrt(pop));
            fprintf(ofid,'color [ORDE]{%i} %6.3f %6.3f %6.3f\n',conformers(c),col);
        end
    end
    fprintf(ofid,'color [DISO]{:}%s darkgrey\n',options.superimposed);
    fprintf(ofid,'color [ORDE]{:}%s darkgrey\n',options.superimposed);
    fprintf(ofid,'view  x\n');
    fprintf(ofid,'detach\n');
    fprintf(ofid,'zoom out\n');
    fprintf(ofid,'copy %s png\n',graphics_x);
    fprintf(ofid,'view  y\n');
    fprintf(ofid,'zoom out\n');
    fprintf(ofid,'copy %s png\n',graphics_y);
    fclose(ofid);
end

function my_map = cluster_color_map

my_map = zeros(101,3);
arg = linspace(0,1,101);
arg = arg';
red_avg = 0;
sig_red = 0.35;
amp_red = 0.75;
green_avg = 0.5;
sig_green = 0.3;
amp_green = 0.7;
blue_avg = 1;
sig_blue = 0.35;
amp_blue = 0.75;
my_map(:,1) = amp_red*exp(-(arg-red_avg).^2/(sqrt(2)*sig_red^2)); 
my_map(:,2) = amp_green*exp(-(arg-green_avg).^2/(sqrt(2)*sig_green^2)); 
my_map(:,3) = amp_blue*exp(-(arg-blue_avg).^2/(sqrt(2)*sig_blue^2)); 

function similarity = get_similarity(D,populations,indices1,indices2)

matches = D(indices1,indices2);
selfmatch1 = D(indices1,indices1);
selfmatch2 = D(indices2,indices2);
var11 = kron(populations(indices1),populations(indices1)').*selfmatch1.^2;
var22 = kron(populations(indices2),populations(indices2)').*selfmatch2.^2;
var12 = kron(populations(indices1),populations(indices2)').*matches.^2;
width1 = sqrt(sum(sum(var11)));
width2 = sqrt(sum(sum(var22)));
similarity = sqrt((width1*width2)/sum(sum(var12)));

