function kappa = discrete_curvature(x, y)
%
% DISCRETE_CURVATURE    Computes discrete curvature of a function given by x
%                       and y value pairs
%
%   kappa = DISCRETE_CURVATURE(x,y)
%
% INPUT
% x             abscissa vector
% y             ordinate vector
%
% OUTPUT
% kappa         discrete curvature
%
% supplied by perplexity.ai, 26.01.2025

    % Compute first derivatives
    dx = gradient(x);
    dy = gradient(y);
    
    % Compute second derivatives
    d2x = gradient(dx);
    d2y = gradient(dy);
    
    % Compute curvature
    kappa = abs(dx.*d2y - dy.*d2x) ./ (dx.^2 + dy.^2).^(3/2);
end