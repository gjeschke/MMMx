function B = reprowvector(VEC,N)
%REPROWVECTOR Replicate a row vector VEC to give an array with N rows.
%   B = reprowvector(VEC,N) creates a large matrix A consisting of an 
%   1-by-N tiling of copies of VEC. 
%
%   Based on built-in REPMAT, but faster for this special case
% 

B = VEC(ones(N, 1), :);


