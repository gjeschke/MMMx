function [dmat,err,res]=metrize_exhaustive(lower_bounds,upper_bounds,pair_indices,fractions,intervals)
% [dmat,err]=metrize_exhaustive(lower_bounds,upper_bounds,pair_indices,fractions,intervals)
%
% Generates a deterministic metric space (distance matrix dmat) for given
% fractions of the range between upper and lower bounds, repeated
% triangular bound smoothing ensures that ranges are as tight as possible
% unspecified restraints are set to the center of the bound range
%
% lower_bounds  lower distance bounds, square matrix [p,p]
% upper_bounds  upper distance bounds, square matrix [p,p]
% pair_indices  index pairs for distances to be set [m,2], 
%               where m = 3 p (p/3-1)/2
% fractions     range fractions where the distances are to be set [m,1]  
% intervals     numbers of interval borders for the restraints, [m,1]
%               number of intervals is larger by one, since the lower and
%               upper bound are not counted and Zaunlattenprinzip applies
%
% dmat          distance matrix [p,p]
% err           error code: 0 no error; 1 erroneous bounds
% res           resolution, maximum interval length
%
% G. Jeschke, 4.10.2017

res = 0;
err=0;

dmat = (lower_bounds+upper_bounds)/2; % initialize distance matrix
[m,~] = size(pair_indices);

for k = 1:m
    [lower_bounds,upper_bounds,err] = triangle_bound_smoothing(lower_bounds,upper_bounds,1);
    if err
        res = [];
        dmat = [];
        return
    end
    p1 = pair_indices(k,1);
    p2 = pair_indices(k,2);
    dmat(p1,p2) = (1-fractions(k))*lower_bounds(p1,p2) + fractions(k)*upper_bounds(p1,p2);
    dmat(p2,p1) = dmat(p1,p2);
    cres = (upper_bounds(p1,p2) - lower_bounds(p1,p2))/(intervals(k)+1);
    if cres > res
        res = cres;
    end
end
