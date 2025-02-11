function [clusters,D] = cluster_transition(entity1,entity2,options)
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
%               .nc             number of clusters, defaults to total
%                               number of conformers over 12
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
%                               also as basis name for the data output file
%               .fname2         basis name for PDB file of ensemble 2, used
%                               also as basis name for the data output file
%               .figname        if present and not empty, the similarity
%                               scaling figure is saved under rhis name
% OUTPUT
% clusters      information on clusters, struct
%               .nc             number of clusters
%               .C1             number of conformers in first ensemble
%               .C2             number of conformers in second ensemble
%               .members        cell(nc,2), members of each cluster
%                               belonging to ensemble 1 and to ensemble 2
%               .type           cluster type, 1 for pure ensemble 1, 2 for
%                               pure ensemble 2, 0 for mixed
%               .similarity     similarity of the two input ensembles
%               .similarities   double matrix (nc,nc) of similarity
%                               measures
%               .subsimilar     double vector (nc,1) of similarities
%                               between subclusters, NaN for pure clusters
%                               or if a cluster contains less then three
%                               conformers of one of the two ensembles
%               .Rg             radius of gyration of cluster
%               .scale          double vector (nc,1) of a one-dimensional
%                               similarity scale between 0 and 1
%               .colorscale     double vector (nc,1) of fractions on a
%                               color scale, ensemble 1 between 0 and 0.2,
%                               mixed ensemble between 0.4 and 0.6,
%                               ensemble 2 between 0.8 and 1
% D     	    dRMSD matrix for the conformers of both ensembles
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

clusters.C1 = length(entity1.populations);
clusters.C2 = length(entity2.populations);

datname = sprintf('transition_%s_%s.csv',options.fname1,options.fname2);
description{1,1} = 'Conformer';
description{2,1} = '%i';
description{1,2} = 'Population';
description{2,2} = '%i';
description{1,3} = 'x';
description{2,3} = '%5.1f';
description{1,4} = 'y';
description{2,4} = '%5.1f';
description{1,5} = 'z';
description{2,5} = '%5.1f';
description{1,6} = 'Uncertainty';
description{2,6} = '%5.1f';
description{1,7} = 'State';
description{2,7} = '%i';
description{1,8} = 'Cluster';
description{2,8} = '%i';
description{1,9} = 'Cluster type';
description{2,9} = '%i';

if ~exist('options','var') || isempty(options) || ~isfield(options,'nc')
    options.nc = round((clusters.C1+clusters.C2)/12);
    % options.nc = 12;
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

if ~isfield(options,'mode')
    options.mode = 'drms';
end

if ~isfield(options,'graphics') || isempty(options.graphics) || isempty(options.graphics{1})
    options.graphics = {};
end

populations = [entity1.populations/sum(entity1.populations);entity2.populations/sum(entity2.populations)];
[D,~,Rg] = pair_drms_matrix(entity1,options.chain1,[],entity2,options.chain2,[]);

m1 = length(entity1.populations);
matches = D(1:m1,m1+1:end);
selfmatch1 = D(1:m1,1:m1);
selfmatch2 = D(m1+1:end,m1+1:end);
var11 = kron(entity1.populations,entity1.populations').*selfmatch1.^2;
var22 = kron(entity2.populations,entity2.populations').*selfmatch2.^2;
var12 = kron(entity1.populations,entity2.populations').*matches.^2;
width1 = sqrt(sum(sum(var11)));
width2 = sqrt(sum(sum(var22)));
clusters.similarity = sqrt((width1*width2)/sum(sum(var12)));

[~,indices] = acs_sorting(D);

[C,~] = size(D);

output_data = zeros(C,9);
output_data(:,1) = 1:C;
output_data(:,2) = populations;
output_data(1:clusters.C1,7) = 1; 
output_data(clusters.C1+1:clusters.C1+clusters.C2,7) = 2; 

sorting = output_data(indices,7);

figure; clf; hold on;
plot(1:C,sorting);

xyz = refined_3D_embedding(D,Rg,populations);
output_data(:,3:5) = xyz; 

D_check = squareform(pdist(xyz));
graces = tinv(0.975,C-1)*sqrt(sum((D_check-D).^2/(C-1)));
output_data(:,6) = graces;

% compute linkage
Z = linkage(D,'complete');
% cluster
assignment = cluster(Z,'maxclust',options.nc);

coor = cmdscale(D);

colors = turbo(4);
colors = 0.75*colors(2:3,:);

[~,score,~,~,explained] = pca(coor);

fh = figure; hold on;
for c = 1:C
    plot(score(c,1),score(c,2),'.','MarkerSize',14,'Color',colors(output_data(c,7),:));
end
axis equal
xlabel('x [Å]');
ylabel('y [Å]');
title(sprintf('PCA 2 analysis (%4.1f%% of variance)',sum(explained(1:2))));
set(gca,'FontSize',12);
[pname,bname,ext] = fileparts(options.figname);
if ~strcmpi(ext,'.tif') && ~strcmpi(ext,'.bmp') && ~strcmpi(ext,'.fig')
    ext = '.pdf';
end
if isfield(options,'figname') && ~isempty(options.figname)
    saveas(fh,fullfile(pname,sprintf('%s_PCA2%s',bname,ext)));
end

fh = figure; hold on;
for c = 1:C
    plot3(score(c,1),score(c,2),score(c,3),'.','MarkerSize',14,'Color',colors(output_data(c,7),:));
end
axis equal
xlabel('x [Å]');
ylabel('y [Å]');
zlabel('z [Å]');
title(sprintf('PCA 3 analysis (%4.1f%% of variance)',sum(explained(1:3))));
set(gca,'FontSize',12);
view(30,30);
if isfield(options,'figname') && ~isempty(options.figname)
    saveas(fh,fullfile(pname,sprintf('%s_PCA3%s',bname,ext)));
end

output_data(:,8) = assignment;

clusters.nc = options.nc;
clusters.members = cell(clusters.nc,2);
clusters.type = zeros(clusters.nc,1);
clusters.Rg = zeros(clusters.nc,1);

h = figure; clf; hold on;
initial_col = [0.2369    0.4437    0.9996];
final_col = [0.9608    0.8902    0.1532];
mixed_col = [0.5 0.5 0.5];
tint = 0.5;
% assign clusters, cluster types 1, 2 and 0 (mixed)
conformers = 1:length(assignment);
for clust = 1:clusters.nc
    members = conformers(assignment == clust);
    clusters.Rg(clust) = sqrt(sum(populations(members).*Rg(members).^2)/sum(populations(members)));
    if max(members) <= clusters.C1 % pure ensemble 1 cluster
        clusters.type(clust) = 1;
        clusters.members{clust,1} = members;
        clusters.members{clust,2} = [];
        for m = 1:length(members)
            pop = populations(members(m))/max(populations);
            conformer_sphere(xyz(members(m),:),graces(members(m)),initial_col,pop);
        end
    end
    if min(members) > clusters.C1 % pure ensemble 2 cluster
        clusters.type(clust) = 2;
        clusters.members{clust,1} = [];
        clusters.members{clust,2} = members - clusters.C1;
        for m = 1:length(members)
            pop = populations(members(m))/max(populations);
            conformer_sphere(xyz(members(m),:),graces(members(m)),final_col,pop);
        end
    end
    if clusters.type(clust) == 0 % mixed cluster
        clusters.members{clust,1} = members(members <= clusters.C1);
        clusters.members{clust,2} = members(members > clusters.C1) - clusters.C1;
        for m = 1:length(members)
            pop = populations(members(m))/max(populations);
            if members(m) <= clusters.C1
                col = mixed_col+tint*initial_col;
            else
                col = mixed_col+tint*final_col;
            end
            conformer_sphere(xyz(members(m),:),graces(members(m)),col,pop);
        end
    end
end
axis equal
my_axes = gca;
view(my_axes,[0,0,1]);
my_axes.CameraUpVector = [1,0,0];
mxyz = min(xyz);
hbar = plot3([mxyz(1)-7.5,mxyz(1)-7.5],[-5,5]+mean(xyz(:,2)),[mxyz(3)-7.5,mxyz(3)-7.5],'k','LineWidth',4);
hlabel = text(mxyz(1)-10,mean(xyz(:,2))+2,mxyz(3)-10,sprintf('10 %c',char(197)),'FontSize',36);
axis off
lighting gouraud
lh = camlight;
pos0 = h.Position;
h.Position = get(0, 'Screensize');
[pname,bname,ext] = fileparts(options.figname);
if ~strcmpi(ext,'.tif') && ~strcmpi(ext,'.bmp') && ~strcmpi(ext,'.fig')
    ext = '.png';
end
if isfield(options,'figname') && ~isempty(options.figname)
    saveas(h,fullfile(pname,sprintf('%s_clustering_x%s',bname,ext)));
end
delete(hbar);
delete(hlabel);
view(my_axes,[0,1,0]);
plot3([mxyz(1)-7.5,mxyz(1)-7.5],[mxyz(2)-7.5,mxyz(2)-7.5],[-5,5]+mean(xyz(:,3)),'k','LineWidth',4);
text(mxyz(1)-10,mxyz(2)-10,mean(xyz(:,3))-2,sprintf('10 %c',char(197)),'FontSize',36);
my_axes.CameraUpVector = [1,0,0];
camlight(lh);
if isfield(options,'figname') && ~isempty(options.figname)
    saveas(h,fullfile(pname,sprintf('%s_clustering_y%s',bname,ext)));
end
h.Position = pos0;

for c = 1:C
    assigned = output_data(c,8);
    output_data(c,9) = clusters.type(assigned);
end

put_csv(datname,output_data,description);

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

if sum(clusters.type) > 0
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
    
    D1 = pair_drms_matrix(entity1,options.chain1);
    [~,indices1] = acs_sorting(D1);
    D2 = pair_drms_matrix(entity2,options.chain2);
    [~,indices2] = acs_sorting(D2);

    dsorting = [indices1,indices2];
    D = D(dsorting,dsorting);
    h = plot_pair_drmsd(D,clusters.C1,clusters.C2);
    [pname,basname,figure_format] = fileparts(options.figname);
    basname = fullfile(pname,basname);
    figname = sprintf('%s_sorted_drms%s',basname,figure_format);
    saveas(h,figname);

    
    clusters.scale = scaled;
    
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
else
    scaled = 0.5*ones(clusters.nc,1);
    clusters.scale = scaled;
    rescaled = scaled;
end

% colors = cluster_color_map;
colors = [initial_col;mixed_col;final_col];
clusters.colorscale = rescaled;


h = figure; clf; hold on;

for clust = 1:clusters.nc
    switch clusters.type(clust)
        case 0
            marker = 'o';
            col = colors(2,:);
        case 1
            marker = 'v';
            col = colors(1,:);
        case 2
            marker = '^' ;
            col = colors(3,:);
    end
    plot(scaled(clust),rescaled(clust),marker,'MarkerSize',10,'Color',col,'MarkerFaceColor',col);
end

set(gca,'FontSize',14);
xlabel('Original similarity scale');
ylabel('Rescaled scale');

if isfield(options,'figname') && ~isempty(options.figname)
    saveas(h,options.figname);
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
    [pname,bname] = fileparts(options.visualize);
    visualize_assignment = fullfile(pname,sprintf('%s_assignment.mmm',bname));
    ofid_assign = fopen(visualize_assignment,'wt');
    fprintf(ofid_assign,'%% MMMx transition visualization script with color by assignment\n');
    fprintf(ofid_assign,'new !\n');
    fprintf(ofid_assign,'pdbload %s\n',options.fname1);
    fprintf(ofid_assign,'pdbload %s\n',options.fname2);
    for clust = 1:clusters.nc
        conformers = clusters.members{clust,1};
        for c = 1:length(conformers)
            switch clusters.type(clust)
                case 0
                    col = mixed_col+tint*initial_col;
                case 1
                    col = initial_col;
                case 2
                    col = final_col;
            end
            maxpop = max(entity1.populations);
            pop = entity1.populations(conformers(c));
            if isempty(options.graphics)
                fprintf(ofid_assign,'show [DISO]{%i} coil %6.3f\n',conformers(c),sqrt(pop/maxpop));
                fprintf(ofid_assign,'color [DISO]{%i} %6.3f %6.3f %6.3f\n',conformers(c),col);
            else
                [ngcmd,~] = size(options.graphics);
                for gcmd = 1:ngcmd
                    fprintf(ofid_assign,'%s [DISO]{%i}(%s)%s %s\n',options.graphics{gcmd,1},conformers(c),options.chain1,options.graphics{gcmd,2},options.graphics{gcmd,3});
                end
                fprintf(ofid_assign,'color [DISO]{%i} %6.3f %6.3f %6.3f\n',conformers(c),col);
                fprintf(ofid_assign,'transparency [DISO]{%i} %5.3f\n',conformers(c),sqrt(pop));
            end
        end
        conformers = clusters.members{clust,2};
        for c = 1:length(conformers)
            switch clusters.type(clust)
                case 0
                    col = mixed_col + tint*final_col;
                case 1
                    col = initial_col;
                case 2
                    col = final_col;
            end
            maxpop = max(entity2.populations);
            pop = entity2.populations(conformers(c));
            if isempty(options.graphics)
                fprintf(ofid_assign,'show [ORDE]{%i} coil %6.3f\n',conformers(c),sqrt(pop/maxpop));
                fprintf(ofid_assign,'color [ORDE]{%i} %6.3f %6.3f %6.3f\n',conformers(c),col);
            else
                [ngcmd,~] = size(options.graphics);
                for gcmd = 1:ngcmd
                    fprintf(ofid_assign,'%s [ORDE]{%i}(%s)%s %s\n',options.graphics{gcmd,1},conformers(c),options.chain2,options.graphics{gcmd,2},options.graphics{gcmd,3});
                end
                fprintf(ofid_assign,'color [ORDE]{%i} %6.3f %6.3f %6.3f\n',conformers(c),col);
                fprintf(ofid_assign,'transparency [ORDE]{%i} %5.3f\n',conformers(c),sqrt(pop));
            end
        end
    end
    fprintf(ofid_assign,'view  x\n');
    fprintf(ofid_assign,'detach\n');
    fprintf(ofid_assign,'zoom out\n');
    fprintf(ofid_assign,'copy %s png\n',fullfile(pname,sprintf('%s_view_x.png',bname)));
    fprintf(ofid_assign,'view  y\n');
    fprintf(ofid_assign,'zoom out\n');
    fprintf(ofid_assign,'copy %s png\n',fullfile(pname,sprintf('%s_view_y.png',bname)));
    fprintf(ofid_assign,'view  -x\n');
    fprintf(ofid_assign,'detach\n');
    fprintf(ofid_assign,'zoom out\n');
    fprintf(ofid_assign,'copy %s png\n',fullfile(pname,sprintf('%s_view_-x.png',bname)));
    fprintf(ofid_assign,'view  -y\n');
    fprintf(ofid_assign,'zoom out\n');
    fprintf(ofid_assign,'copy %s png\n',fullfile(pname,sprintf('%s_view_-y.png',bname)));
    fprintf(ofid_assign,'view  z\n');
    fprintf(ofid_assign,'detach\n');
    fprintf(ofid_assign,'zoom out\n');
    fprintf(ofid_assign,'copy %s png\n',fullfile(pname,sprintf('%s_view_z.png',bname)));
    fprintf(ofid_assign,'view  -z\n');
    fprintf(ofid_assign,'zoom out\n');
    fprintf(ofid_assign,'copy %s png\n',fullfile(pname,sprintf('%s_view_-z.png',bname)));
    fclose(ofid_assign);
end

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

function conformer_sphere(coor,grace,col,pop)
[x,y,z,t] = point2trisphere(coor,grace/2,2);
obj = trisurf(t,x,y,z);
set(obj, 'FaceColor', col, 'EdgeColor', 'none','FaceLighting','gouraud','Clipping','off');
set(obj, 'CDataMapping','direct','AlphaDataMapping','none','FaceAlpha',pop);

function [x,y,z,t] = point2trisphere(coor,radius,n)
% Creates a triangulated sphere centered at coordinates coor with radius r
% function is just a wrapper for BuildSphere by Luigi Giaccari
%
% coor      [x0,y0,z0] coordinate vector of the center of the sphere
% radius    radius of the sphere
% n         optional quality number, between 0 and inf, be careful, number
%           of points grows dramatically with increasing n, defaults to 1 
%
% x,y,z     coordinates of triangle vertices
% t         face matrix for calling trisurf, use as
%           trisurf(t,p(:,1),p(:,2),p(:,3))
%
% G. Jeschke, 2009

if ~exist('n','var') || isempty(n)
    n = 2;
end

[p,t] = BuildSphere(n);

p=radius*p;
offset=ones(size(p(:,1)));
x=p(:,1)+coor(1)*offset;
y=p(:,2)+coor(2)*offset;
z=p(:,3)+coor(3)*offset;
function [p,t]=BuildSphere(n)
%BUILDSPHERE: Returns the triangulated model of a sphere using the
%icosaedron subdivision method.
%
%INPUT: 
%n (integer number) indicates the number of subdivisions,
% it can assumes values between 0-inf. The greater n the better will look
%the surface but the more time will be spent in surface costruction and
%more triangles will be put in the output model.
%
%OUTPUT:
%In p (n x 3) and t(m x 3) we can find points and triangles indexes
%of the model. The sphere is supposed to be of unit radius and centered in
%(0,0,0). To obtain spheres centered in different location, or with
%different radius, is just necessary a traslation and a scaling
%trasformation. These operation are not perfomed by this code beacuse it is
%extrimely convinient, in order of time perfomances, to do this operation
%out of the function avoiding to call the costruction step each time.
%
%NOTE: 
%This function is more efficient than the matlab command sphere in
%terms of dimension fo the model/ accuracy of recostruction. This due to
%well traingulated model that requires a minor number of patches for the
%same geometrical recostruction accuracy. Possible improvement are possible
%in time execution and model subdivision flexibilty.
%
%EXAMPLE:
%
% n=5;
%
% [p,t]=BuildSphere(n);
%
% figure(1) axis equal hold on trisurf(t,p(:,1),p(:,2),p(:,3)); axis vis3d
% view(3)
%
%
%CONTACTS:
%
%Author: Giaccari Luigi Created:25/04/2009
%
%For info/bugs/questions/suggestions : giaccariluigi@msn.com
%
%Visit: http://giaccariluigi.altervista.org/blog/

%errors check
if nargin>2 || nargin<1
    error('Wrong number of inputs: possible choice are two or one');
end

if nargout>2 
    error('Maximum two outputs supported');
end

if not(isscalar(n))
    error('input n must be a scalar');
end

if round(n)~=n
    error('input n must be a integer value (float or integer type)');
end

%Coordinates have benn taken from Jon Leech files

% Twelve vertices of icosahedron on unit sphere
tau = 0.8506508084; % t=(1+sqrt(5))/2, tau=t/sqrt(1+t^2)
one = 0.5257311121; % one=1/sqrt(1+t^2) , unit sphere

p( 1,:) = [  tau,  one,    0 ]; % ZA
p( 2,:) = [ -tau,  one,    0 ]; % ZB
p( 3,:) = [ -tau, -one,    0 ]; % ZC
p( 4,:) = [  tau, -one,    0 ]; % ZD
p( 5,:) = [  one,   0 ,  tau ]; % YA
p( 6,:) = [  one,   0 , -tau ]; % YB
p( 7,:) = [ -one,   0 , -tau ]; % YC
p( 8,:) = [ -one,   0 ,  tau ]; % YD
p( 9,:) = [   0 ,  tau,  one ]; % XA
p(10,:) = [   0 , -tau,  one ]; % XB
p(11,:) = [   0 , -tau, -one ]; % XC
p(12,:) = [   0 ,  tau, -one ]; % XD

% Structure for unit icosahedron
t = [  5,  8,  9 ;
    5, 10,  8 ;
    6, 12,  7 ;
    6,  7, 11 ;
    1,  4,  5 ;
    1,  6,  4 ;
    3,  2,  8 ;
    3,  7,  2 ;
    9, 12,  1 ;
    9,  2, 12 ;
    10,  4, 11 ;
    10, 11,  3 ;
    9,  1,  5 ;
    12,  6,  1 ;
    5,  4, 10 ;
    6, 11,  4 ;
    8,  2,  9 ;
    7, 12,  2 ;
    8, 10,  3 ;
    7,  3, 11 ];

nt=20;%we now have 20 triangles
np=12;
%get the final number of points
totp=np;
for i=1:n
    totp=(totp*4)-6;
end
p=[p;zeros(totp-12,3)];%points array

for i=1:n
    
    told=t;
    t=zeros(nt*4,3);%each iteration the number of traingles increase
    peMap=sparse([],[],[],np,np,np);%Point edge sparse matrix
    %(in future try with a vector structure)
    
    ct=1;%counter triangles
    for j=1:nt%loop trough all old triangles
        

        %sort t index (in future use binar logic instead of sort
        %function)

        p1=told(j,1);p2=told(j,2);p3=told(j,3);

        %make temp scalar
        x1=p(p1,1);x2=p(p2,1);x3=p(p3,1);
        y1=p(p1,2);y2=p(p2,2);y3=p(p3,2);
        z1=p(p1,3);z2=p(p2,3);z3=p(p3,3);
        
        
        %first edge
        
        %this line simplify acces to sparse matrix preserving triangles orientation
        if p1<p2;p1m=p1;
            p2m=p2;
        else 
            p2m=p1;
            p1m=p2;
        end
        
        id=peMap(p1m,p2m);
        if id>0%point exist
            p4=id;
        else
            np=np+1;
            p4=np;
            peMap(p1m,p2m)=np;%#ok<SPRIX> %update map
            %get a new point, normalize to one
            x=(x1+x2)/2;
            y=(y1+y2)/2;
            z=(z1+z2)/2;
            l=sqrt(x*x+y*y+z*z);
            x=x/l;y=y/l;z=z/l;
            p(np,1)=x; p(np,2)=y; p(np,3)=z;
        end
        
        
        %second edge
        
        %this line simplify acces to sparse matrix preserving triangles
        %orientation
        if p2<p3
            p2m=p2;
            p3m=p3;
        else
            p2m=p3;
            p3m=p2;
        end
        
        id=peMap(p2m,p3m);
        if id>0%point exist
            p5=id;
        else
            np=np+1;
            p5=np;
            peMap(p2m,p3m)=np;%#ok<SPRIX> %update map
            %get a new point, normalize to one
            x=(x2+x3)/2;
            y=(y2+y3)/2;
            z=(z2+z3)/2;
            l=sqrt(x*x+y*y+z*z);
            x=x/l;y=y/l;z=z/l;
            p(np,1)=x; p(np,2)=y; p(np,3)=z;
        end
        
        
        %third edge
        
        %this line simplify acces to sparse matrix preserving triangles
        %orientation
        if p1<p3
            p1m=p1;
            p3m=p3;
        else
            p3m=p1;
            p1m=p3;
        end
        
        id=peMap(p1m,p3m);
        if id>0%point exist
            p6=id;
        else
            np=np+1;
            p6=np;
            peMap(p1m,p3m)=np;%#ok<SPRIX> %update map
            
            %get a new point, normalize to one
            x=(x1+x3)/2;
            y=(y1+y3)/2;
            z=(z1+z3)/2;
            l=sqrt(x*x+y*y+z*z);
            x=x/l;y=y/l;z=z/l;
            p(np,1)=x; p(np,2)=y; p(np,3)=z;
            
        end
        
        %allocate new triangles
        %               refine indexing
        %                        p1
        %                       /\
        %                      /t1\
        %                   p6/____\p4
        %                    /\    /\
        %                   /t4\t2/t3\
        %                  /____\/____\
        %                 p3	p5     p2
        
        t(ct,1)=p1; t(ct,2)=p4; t(ct,3)=p6;
        ct=ct+1;
        t(ct,1)=p4; t(ct,2)=p5; t(ct,3)=p6;
        ct=ct+1;
        t(ct,1)=p4; t(ct,2)=p2; t(ct,3)=p5;
        ct=ct+1;
        t(ct,1)=p6; t(ct,2)=p5; t(ct,3)=p3;
        ct=ct+1;
        
    end
    
    nt=ct-1;
    
end

function h = plot_pair_drmsd(D,C1,C2)

h = figure; hold on

plot([0,0],[1,C1],'Color',[0.75,0,0],'LineWidth',2);
plot([0,0],[C1+1,C1+C2],'Color',[0,0,0.8],'LineWidth',2);
plot([1,C1],[0,0],'Color',[0.75,0,0],'LineWidth',2);
plot([C1+1,C1+C2],[0,0],'Color',[0,0,0.8],'LineWidth',2);
image(D,'CDataMapping','scaled');
curr_axis = gca;
curr_axis.YDir = 'normal';
colorbar;
axis tight
xlabel('Conformer number');
ylabel('Conformer number');
title('Sorted distance root mean square deviation');
axis equal
