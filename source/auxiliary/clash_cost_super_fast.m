function cost = clash_cost_super_fast(ar,kpa,br,kpb,a,b,ct,ft,tol)
% cost = clash_cost_super_fast(ar,kpa,br,kpb,a,b,ct,ft,tol)
%
% Computes the clash cost of two convex hulls defined as their overlap
% volume, points are identified that are common two both hulls, allowing
% for some tolerance (default of 1.5 Å is a good choice for macromolecules)
% the volume of the convex hull for these points is the clash cost
% uses inhull by John D'Errico
% use convhulln to generate the input from atom coordinates
%
% ar    Cartesian coordinates of vertices of the first hull
% kpa   indices defining triangular faces of the first hull
% br    Cartesian coordinates of vertices of the second hull
% kpb   indices defining triangular faces of the second hull
% a     coordinates of all heavy atoms of first molecule
% b     coordinates of all heavy atoms of second molecule
% ct    clash threshold for calling clash_cost if required
% ft    fine threshold, cost below which the more expensive computation is
%       done, defaults to 200
% tol   tolerance for inclusion of a point in a convex hull, defaults to
%       1.5
%
% cost  overlap volume
%
% G. Jeschke, 10.8.2018

cost = 0;

if ~exist('ft','var') || isempty(ft)
    ft = 200; % if clash cost is non-zero, but drops below ft, the slow 
              % fine-grained clash_cost function is called
end

if ~exist('tol','var') || isempty(tol)
    tol = 1.5;
end

inb = inhull(ar,br,kpb,tol);
p1 = ar(inb,:);
ina = inhull(br,ar,kpa,tol);
p2 = br(ina,:);
p = [p1;p2];

[np,~] = size(p);
if np > 3
    try
        [~,cost] = convhulln(double(p));
    catch
        cost = 1e6;
    end
end

if cost > eps && cost < ft
    cost = clash_cost(a,b,ct);
end