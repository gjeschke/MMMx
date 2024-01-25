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
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

clusters.C1 = length(entity1.populations);
clusters.C2 = length(entity2.populations);

if ~exist('options','var') || isempty(options) || ~isfield(options,'nc')
    % options.nc = min([round(clusters.C1/10),round(clusters.C2/10)]);
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

if ~isfield(options,'mode')
    options.mode = 'drms';
end

if ~isfield(options,'graphics') || isempty(options.graphics) || isempty(options.graphics{1})
    options.graphics = {};
end

populations = [entity1.populations;entity2.populations];
[D,~,Rg] = pair_drms_matrix(entity1,options.chain1,[],entity2,options.chain2,[]);

[C,~] = size(D);
% severity = zeros(1,3*C*(C-1)*(C-2));
% triples = zeros(3*C*(C-1)*(C-2),3);
% violations = 0;
% for c1 = 1:C-2
%     for c2 = c1+1:C-1
%         for c3 = c2+1:C
%             d12 = D(c1,c2);
%             d23 = D(c2,c3);
%             d13 = D(c1,c3);
%             if d12 > d13 + d23
%                 violations = violations + 1;
%                 severity(violations) = d12 - (d13+d23);
%                 triples(violations,:) = [c1,c2,c3];
%             end
%             if d13 > d12 + d23
%                 violations = violations + 1;
%                 severity(violations) = d13 - (d12+d23);
%                 triples(violations,:) = [c1,c2,c3];
%             end
%             if d23 > d12 + d13
%                 violations = violations + 1;
%                 severity(violations) = d23 - (d12+d13);
%                 triples(violations,:) = [c1,c2,c3];
%             end
%         end
%     end
% end
% 
% severity = severity(1:violations);
% triples = triples(1:violations,:);
% [maxsev,index] = max(severity);

% c = triples(index,:);
% fprintf(1,'%i violations of the triangle inequality were found\n',violations);
% if violations > 0
%     fprintf(1,'Maximum severity is %4.1f %c\n',maxsev,char(197));
%     fprintf(1,'Maximum severity for triple (%i,%i,%i)\n',c);
%     fprintf(1,'The conformer distances are (%4.1f,%4.1f,%4.1f) %c\n',pdm(c(1),c(2)),pdm(c(1),c(3)),pdm(c(2),c(3)),char(197));
% end

coor0 = dmat2coor(D);
if isempty(coor0)
    fprintf(1,'DRMSD matrix could not be embedded into 3D space\n');
else
    [coor1,grace,all_violations] = adaptive_embedding_refiner(coor0,D,1);
    fprintf(1,'Embedding was refined with maximum violation %4.2f %c\n',grace,char(197));
    figure(333); clf; hold on;
    plot(all_violations(all_violations > 0),'.','MarkerSize',10,'Color',[0.75,0,0]);
end

% get the inertia tensor in the original frame
centred = true;
inertia = get_inertia_tensor(coor1,centred);
% diagonalize the inertia tensor
[evec,ID] = eig(inertia);
% sort eigenvalues in ascending order
[~,indices] = sort(diag(ID));
% sort eigenvectors accordingly
evecs = evec(:,indices);
% determine the centre of gravity of the conformer set
xyzc = mean(coor1);
% put orgin at centre of gravity of the conformer set
xyz = coor1 - repmat(xyzc,C,1);
% transform into inertia frame
xyz = xyz*evecs;

[Rg_sorted,Rg_sorting] = sort(Rg);
xsorting = xyz(Rg_sorting,1);
invertx = sum((1:C)'.*xsorting) < 0;
ysorting = xyz(Rg_sorting,2);
inverty = sum((1:C)'.*ysorting) < 0;
zsorting = xyz(Rg_sorting,3);

figure(111); clf; hold on
plot(Rg_sorted,'k.');

figure(222); clf; hold on;
plot(xsorting,'.','Color',[0.75,0,0]);
plot(ysorting,'.','Color',[0,0.6,0]);
plot(zsorting,'.','Color',[0,0,0.8]);

% invert coordinates so that lowest Rg terminus has lower x and y coordinates
if invertx % the x coordinate must be inverted
    xyz(:,1) = -xyz(:,1);
    % this requires that either the y or the z coordinate is also inverted
    % (keep chirality)
    if inverty % in this case, the y coordinate should be inverted
        xyz(:,2) = -xyz(:,2);
    else % otherwise, the z coordinate is inverted
        xyz(:,3) = -xyz(:,3);
    end
else % the x coordinate was not inverted
    % if the y coordinate must be inverted, the z coordinate must be
    % inverted, too, to keep the frame right-handed
    if inverty % in this case, the x coordinate should be inverted
        xyz(:,2) = -xyz(:,2);
        xyz(:,3) = -xyz(:,3);
    end
end

% compute linkage
Z = linkage(D,'complete');
% cluster
assignment = cluster(Z,'maxclust',options.nc);
clusters.nc = options.nc;
clusters.members = cell(clusters.nc,2);
clusters.type = zeros(clusters.nc,1);
clusters.Rg = zeros(clusters.nc,1);

h = figure; clf; hold on;
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
            col = [0.75,0,0];
            conformer_sphere(xyz(members(m),:),grace,col);
        end
    end
    if min(members) > clusters.C1 % pure ensemble 2 cluster
        clusters.type(clust) = 2;
        clusters.members{clust,1} = [];
        clusters.members{clust,2} = members - clusters.C1;
        for m = 1:length(members)
            col = [0,0,0.75];
            conformer_sphere(xyz(members(m),:),grace,col);
        end
    end
    if clusters.type(clust) == 0 % mixed cluster
        clusters.members{clust,1} = members(members <= clusters.C1);
        clusters.members{clust,2} = members(members > clusters.C1) - clusters.C1;
        for m = 1:length(members)
            if members(m) <= clusters.C1
                col = [0.42,0.56,0.14];
            else
                col = [0.18,0.55,0.34];
            end
            conformer_sphere(xyz(members(m),:),grace,col);
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
    
    [~,csorting] = sort(scaled);
    all_indices1 = zeros(1,clusters.C1);
    poi1 = 0;
    all_indices2 = zeros(1,clusters.C2);
    poi2 = 0;
    for kc = 1:length(csorting)
        clust = csorting(kc);
        indices1 = clusters.members{clust,1};
        all_indices1(poi1+1:poi1+length(indices1)) = indices1;
        poi1 = poi1+length(indices1);
        indices2 = clusters.members{clust,2} + clusters.C1;
        all_indices2(poi2+1:poi2+length(indices2)) = indices2;
        poi2 = poi2+length(indices2);
    end
    dsorting = [all_indices1,all_indices2];
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

colors = cluster_color_map;
clusters.colorscale = rescaled;


h = figure; clf; hold on;

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

set(gca,'FontSize',14);
xlabel('original similarity scale');
ylabel('rescaled values for color code');

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
    visualize_cs = fullfile(pname,sprintf('%s_conformational_selection.mmm',bname));
    ofid_cs = fopen(visualize_cs,'wt');
    visualize_dp = fullfile(pname,sprintf('%s_depopulated.mmm',bname));
    ofid_dp = fopen(visualize_dp,'wt');
    visualize_if = fullfile(pname,sprintf('%s_induced_fit.mmm',bname));
    ofid_if = fopen(visualize_if,'wt');
    fprintf(ofid_cs,'%% MMMx transition visualization script for conformational selection\n');
    fprintf(ofid_cs,'new !\n');
    fprintf(ofid_cs,'pdbload %s\n',options.fname1);
    fprintf(ofid_cs,'pdbload %s\n',options.fname2);
    fprintf(ofid_dp,'%% MMMx transition visualization script for depopulated conformers\n');
    fprintf(ofid_dp,'new !\n');
    fprintf(ofid_dp,'pdbload %s\n',options.fname1);
    fprintf(ofid_dp,'pdbload %s\n',options.fname2);
    fprintf(ofid_if,'%% MMMx transition visualization script for induced fit\n');
    fprintf(ofid_if,'new !\n');
    fprintf(ofid_if,'pdbload %s\n',options.fname1);
    fprintf(ofid_if,'pdbload %s\n',options.fname2);
    for clust = 1:clusters.nc
        switch clusters.type(clust)
            case 0
                ofid = ofid_cs;
                col = [0,0.6,0];
            case 1
                ofid = ofid_dp;
                col = [0.75,0,0];
            case 2
                ofid = ofid_if;
                col = [0,0,0.75];
        end
%         color_index = 1 + round(100*rescaled(clust));
%         col = colors(color_index,:);
        conformers = clusters.members{clust,1};
        for c = 1:length(conformers)
            pop = entity1.populations(conformers(c))/max(entity1.populations);
            if isempty(options.graphics)
                fprintf(ofid,'show [DISO]{%i} coil %6.3f\n',conformers(c),sqrt(pop));
                fprintf(ofid,'color [DISO]{%i} %6.3f %6.3f %6.3f\n',conformers(c),col);
            else
                [ngcmd,~] = size(options.graphics);
                for gcmd = 1:ngcmd
                    fprintf(ofid,'%s [DISO]{%i}(%s)%s %s\n',options.graphics{gcmd,1},conformers(c),options.chain1,options.graphics{gcmd,2},options.graphics{gcmd,3});
                end
                fprintf(ofid,'transparency [DISO]{%i} %5.3f\n',conformers(c),pop);
            end
        end
        conformers = clusters.members{clust,2};
        for c = 1:length(conformers)
            pop = entity2.populations(conformers(c))/max(entity2.populations);
            if isempty(options.graphics)
                fprintf(ofid,'show [ORDE]{%i} coil %6.3f\n',conformers(c),sqrt(pop));
                fprintf(ofid,'color [ORDE]{%i} %6.3f %6.3f %6.3f\n',conformers(c),col);
            else
                [ngcmd,~] = size(options.graphics);
                for gcmd = 1:ngcmd
                    fprintf(ofid,'%s [ORDE]{%i}(%s)%s %s\n',options.graphics{gcmd,1},conformers(c),options.chain2,options.graphics{gcmd,2},options.graphics{gcmd,3});
                end
                fprintf(ofid,'transparency [ORDE]{%i} %5.3f\n',conformers(c),pop);
            end
        end
    end
    if isempty(options.graphics)
        fprintf(ofid_cs,'color [DISO]{:}%s darkgrey\n',options.superimposed);
        fprintf(ofid_cs,'color [ORDE]{:}%s darkgrey\n',options.superimposed);
        fprintf(ofid_dp,'color [DISO]{:}%s darkgrey\n',options.superimposed);
        fprintf(ofid_dp,'color [ORDE]{:}%s darkgrey\n',options.superimposed);
        fprintf(ofid_if,'color [DISO]{:}%s darkgrey\n',options.superimposed);
        fprintf(ofid_if,'color [ORDE]{:}%s darkgrey\n',options.superimposed);
    end
    types = {'conformational_selection' 'depopulated' 'induced_fit'};
    ofids = [ofid_cs,ofid_dp,ofid_if];
    for tp = 1:3
        ofid = ofids(tp);
        my_type = types{tp};
        fprintf(ofid,'view  x\n');
        fprintf(ofid,'detach\n');
        fprintf(ofid,'zoom out\n');
        fprintf(ofid,'copy %s png\n',fullfile(pname,sprintf('%s_%s_view_x.png',bname,my_type)));
        fprintf(ofid,'view  y\n');
        fprintf(ofid,'zoom out\n');
        fprintf(ofid,'copy %s png\n',fullfile(pname,sprintf('%s_%s_view_y.png',bname,my_type)));
        fprintf(ofid,'view  -x\n');
        fprintf(ofid,'detach\n');
        fprintf(ofid,'zoom out\n');
        fprintf(ofid,'copy %s png\n',fullfile(pname,sprintf('%s_%s_view_-x.png',bname,my_type)));
        fprintf(ofid,'view  -y\n');
        fprintf(ofid,'zoom out\n');
        fprintf(ofid,'copy %s png\n',fullfile(pname,sprintf('%s_%s_view_-y.png',bname,my_type)));
        fprintf(ofid,'view  z\n');
        fprintf(ofid,'detach\n');
        fprintf(ofid,'zoom out\n');
        fprintf(ofid,'copy %s png\n',fullfile(pname,sprintf('%s_%s_view_z.png',bname,my_type)));
        fprintf(ofid,'view  -z\n');
        fprintf(ofid,'zoom out\n');
        fprintf(ofid,'copy %s png\n',fullfile(pname,sprintf('%s_%s_view_-z.png',bname,my_type)));
    end
    fclose(ofid_cs);
    fclose(ofid_dp);
    fclose(ofid_if);
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

function conformer_sphere(coor,grace,col)
[x,y,z,t] = point2trisphere(coor,grace/2,2);
obj = trisurf(t,x,y,z);
set(obj, 'FaceColor', col, 'EdgeColor', 'none','FaceLighting','gouraud','Clipping','off');
set(obj, 'CDataMapping','direct','AlphaDataMapping','none');

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
