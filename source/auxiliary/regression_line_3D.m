function [p0,v] = regression_line_3D(xyz)
% function [p0,v] = regression_line_3D(xyz)
%
% linear regression line through points in 3D space
%
% xyz   Cartesian coordinates of the points nx3 array
%
% p0    a point on the line, midpoint of point cloud
% v     vector along the linear regression line, length approximates
%       extension of the point cloud along this direction
%
% G. Jeschke, 2009-2021

% make sure that the line intersects the xy plane (subtract center
% coordinate)
xyzc = mean(xyz);
[n,~] = size(xyz);
xyz = xyz - repmat(xyzc,n,1);

dx = max(xyz(:,1)) - min(xyz(:,1));
dy = max(xyz(:,2)) - min(xyz(:,2));
dz = max(xyz(:,3)) - min(xyz(:,3));

[~,indy] = sort([dx,dy,dz]);

xyz1 = [xyz(:,indy(1)) xyz(:,indy(2)) xyz(:,indy(3))];

% linear 3D regression, see http://en.wikipedia.org/wiki/User:Vossman/3D_Line_Regression

sx = sum(xyz1(:,1));
sy = sum(xyz1(:,2));
sz = sum(xyz1(:,3));
sz2 = sum(xyz1(:,3).^2);
sxz = sum(xyz1(:,1).*xyz1(:,3));
syz = sum(xyz1(:,2).*xyz1(:,3));

x0 = (sxz*sz-sx*sz2)/(sz^2-n*sz2);
y0 = (syz*sz-sy*sz2)/(sz^2-n*sz2);
vx = (n*sxz-sx*sz)/(n*sz2-sz^2);
vy = (n*syz-sy*sz)/(n*sz2-sz^2);

p0 = zeros(1,3);
p0(indy(1)) = x0;
p0(indy(2)) = y0;
p0 = p0 + xyzc;
v = ones(1,3);
v(indy(1)) = vx;
v(indy(2)) = vy;
v = v/norm(v);

dd = sqrt(dx^2+dy^2+dz^2);
v = dd*v;
