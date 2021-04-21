function [rmsd, coor2b, transmat] = superimpose_3points(coor1,coor2)
%
% [rmsd, coor2b, transmat] = superimpose_3points(coor1,coor2)
%
%  Standard-frame based superposition of three points in space
%  unlike least-squares fitting, this is unambiguous by ensuring that the
%  standard frame is a right-handed frame, the sequence of the points
%  resolves ambiguity
%  the second point (frame origin) is matched exactly
% 
% equal weights wi=1/m are assumed for all atomic positions
%
% coor1   master coordinate set, must be an (m x 3) array
% coor2   coordinate set to be fitted to coor1 by translation and rotation,
%         must be an (m x 3) array with the same m, otherwise rms=-1e6 is
%         output
% weights vector of length m with weights for the atomic positions, can be
%         column or row vector, defaults to uniform weights
%
% rms     root mean square deviation of both coordinate sets
% coor2b  coordinate set 2 transformed to give the best fit with coordinate
%         set 1
% transmat   affine transformation matrix for coordinate set 2
%
% transmat is supplied so that after a fit of a substructure other parts of
% structure 2 can be transformed with transform_structure (or simple
% multiplication of a coordinate array with transmat)
%
% (c) G. Jeschke, 12.1.2017


[Rp1,orig1] = get_trafo(coor1);
[Rp2,orig2] = get_trafo(coor2);

coor2a = coor2 - repmat(orig2,3,1);
coor2a = coor2a*Rp2';
coor2b = coor2a*Rp1 + repmat(orig1,3,1);
rmsd = sqrt(sum(sum((coor2b-coor1).^2))/3);
transmat1a = zeros(4);
transmat1a(1:3,1:3) = diag([1,1,1]);
transmat1a(1:3,4) = -orig2';
transmat1a(4,4) = 1;
% coor2ax = [coor2 ones(3,1)]*transmat1a';
transmat1b = zeros(4);
transmat1b(1:3,1:3) = Rp2;
transmat1b(4,4) = 1;
% coor2ax = coor2ax*transmat1b';
transmat2 = zeros(4);
transmat2(1:3,1:3) = Rp1';
transmat2(1:3,4) = orig1';
transmat2(4,4) = 1;
% coor2bx = coor2ax*transmat2';
transmat = transmat2*transmat1b*transmat1a;
% coor2c = [coor2 ones(3,1)]*transmat';

function [Rp,orig] = get_trafo(coor)

orig = coor(2,:);
coor = coor - repmat(orig,3,1);
x = coor(1,:)-coor(2,:); 
x = x/norm(x);    % unit vector along x
yp = coor(3,:)-coor(2,:); 
yp = yp/norm(yp);
z = cross_rowvec(x,yp); % z axis is perpendicular on xy plane
z = z/norm(z);
y = cross_rowvec(z,x); % real (corrected) y axis
Rp = [x;y;z];