function [mean_val, median_val] = compute_mean_median(x, px, level)
% Compute mean and median from probability density function
% if input 'level' is provided and not 0.5, median_val is replaced by the
% value xl where probability p = level for the integratl from min(x) to xl
% x: vector of axis values
% px: vector of probability density values
% level: level for 'median' value, defaults to 0.5

px = px/sum(px);

mean_val = sum(x .* px);
cdf = cumsum(px);
[~,idx] = min(abs(cdf-level));
range = idx - 1 : idx+1;
dx = (x(idx+1)-x(idx-1))/2;
x0 = x(range);
xf = x(idx-1):dx/10:x(idx+1);
cdff = interp1(x0,cdf(range),xf,"linear",'extrap');
[~,idx] = min(abs(cdff-level));
median_val = xf(idx);

end