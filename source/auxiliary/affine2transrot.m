function [trans,euler] = affine2transrot(transmat)  
% Extracts the translation vector and Euler angles from an affine
% transformation matrix, assuming that the rotation is performed first and
% the translation second
%
% see: https://gamedev.stackexchange.com/questions/50963/how-to-extract-euler-angles-from-transformation-matrix
%
% G. Jeschke, 4.10.2017

euler = zeros(1,3);
euler(1)  = atan2(-transmat(2,3), transmat(3,3));
cx = cos(euler(1));
sx = sin(euler(1));
euler(2) = asin(transmat(1,3));
sz = cx*transmat(2,1) + sx*transmat(3,1);
cz = cx*transmat(2,2) + sx*transmat(3,2);
euler(3) = atan2(sz, cz);

trans = transmat(1:3,4)';
