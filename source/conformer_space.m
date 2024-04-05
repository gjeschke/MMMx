function [conf_space,assignment] = conformer_space(entity,options)

conf_space = [];

if ~exist('options','var') || isempty(options) || ~isfield(options,'mode')
    options.mode = 'drms';
end

switch options.mode
    case 'drms'
        [D,~,Rg] = pair_drms_matrix(entity);
    case 'rmsd'
        [D,~,~,Rg] = pair_rmsd_matrix(entity);
    case 'rmsd_oriented'
        [D,~,~,Rg] = pair_rmsd_matrix_oriented(entity);
    otherwise
        return
end

[~,~,assignment,cluster_sizes] = cluster_sorting(D,entity.populations);

col = turbo(length(cluster_sizes));

[C,~] = size(D);

coor0 = dmat2coor(D);
if isempty(coor0)
    fprintf(1,'DRMSD matrix could not be embedded into 3D space\n');
else
    [coor1,grace] = adaptive_embedding_refiner(coor0,D,1);
    fprintf(1,'Embedding was refined with maximum violation %4.2f %c\n',grace,char(197));
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

[~,Rg_sorting] = sort(Rg);
xsorting = xyz(Rg_sorting,1);
invertx = sum((1:C)'.*xsorting) < 0;
ysorting = xyz(Rg_sorting,2);
inverty = sum((1:C)'.*ysorting) < 0;

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

conf_space = xyz;

h = figure; clf; hold on;
% assign clusters, cluster types 1, 2 and 0 (mixed)
for c = 1:C
    conformer_sphere(xyz(c,:),grace,col(assignment(c),:));
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
    saveas(h,fullfile(pname,sprintf('%s_conformer_space_x%s',bname,ext)));
end
delete(hbar);
delete(hlabel);
view(my_axes,[0,1,0]);
plot3([mxyz(1)-7.5,mxyz(1)-7.5],[mxyz(2)-7.5,mxyz(2)-7.5],[-5,5]+mean(xyz(:,3)),'k','LineWidth',4);
text(mxyz(1)-10,mxyz(2)-10,mean(xyz(:,3))-2,sprintf('10 %c',char(197)),'FontSize',36);
my_axes.CameraUpVector = [1,0,0];
camlight(lh);
if isfield(options,'figname') && ~isempty(options.figname)
    saveas(h,fullfile(pname,sprintf('%s_conformer_space_y%s',bname,ext)));
end
h.Position = pos0;

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
