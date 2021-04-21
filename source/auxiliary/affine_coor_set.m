function xyz = affine_coor_set(xyz,transmat)
% function xyz = affine_coor_set(xyz,transmat)
%
% affine coordinate transformation of an m*3 set of Cartesian point
% coordinates using the 4*4 affine transformation matrix transmat or the
% set of affine transforamtion matrices that are applied consecutively
%
% G. Jeschke, 2010

if iscell(transmat)
    matrix=eye(4);
    num=length(transmat);
    for k=1:num
        matrix=transmat{k}*matrix;
    end
else
    matrix=transmat;
end

[m,~]=size(xyz);
xyz=[xyz ones(m,1)];
xyz=matrix*xyz';
xyz=xyz';
xyz=xyz(:,1:3);
