function cost = clash_cost(a,b,threshold,tolerance)
% cost = clash_cost(a,b,threshold,tolerance)
%
% Returns cost function value for clashes between atoms with coordinates 
% in arrays a and b
%
% atoms closer than the input threshold are considered as clashing
%
% constant tolerance (optional) defines how steeply cost increases with
% violation (it is a standard deviation for a Gaussian tolerance function)
% in Ångstroem, defaults to 0.1 Ångstroem
%
% cost = clash_cost(a,b,threshold)
%
% G. Jeschke, 9.11.2017

if threshold >1e4
    cost = 0;
    return
end
if ~exist('tolerance','var')
    tolerance = single(0.1);
end


[m1,~] = size(a);
[m2,~] = size(b);
a2 = repmat(sum(a.^2,2),1,m2);
b2 = repmat(sum(b.^2,2),1,m1).';
pair_dist = sqrt(abs(a2 + b2 - 2*a*b.'));
clashing = pair_dist(pair_dist <= threshold);
cost = sum((threshold-clashing).^2)/tolerance^2;

% fprintf(1,'Cost: %12.1f\n',cost);