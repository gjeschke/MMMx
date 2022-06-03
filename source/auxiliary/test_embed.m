function [lower_dev,upper_dev]=test_embed(coor,lower_bounds,upper_bounds)
%[lower_dev,upper_dev]=test_embed(coor,lower_bounds,upper_bounds)
%
% Test for constraint violations of an embedded structure with coordinates
% coor against lower and upper bounds for all distances
%
% lower_dev   matrix of violations of lower bounds
% upper_dev   matrix of violations of upper bounds
%
% (c) G. Jeschke, 2008

dmat=coor2dmat(coor);

lower_dev=lower_bounds-dmat;
lower_dev=lower_dev.*(lower_dev>0);

upper_dev=dmat-upper_bounds;
upper_dev=upper_dev.*(upper_dev>0);

