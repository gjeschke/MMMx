function [points,mindist] = fibonacci_sphere(n)
% [points,mindist] = fibonacci_sphere(n)
%
% Provides a list of n points on a unit sphere in Cartesian coordinate form
% that are approximately uniformly distributed (using a Fibonacci spiral),
% the minimum pair distance mindist can be requested
%
% G. Jeschke, 10.9.2019

points = zeros(n,3);
offset = 2/n;
increment = pi * (3 - sqrt(5));

for k = 1:n
    z = (((k-1) * offset) - 1) + (offset / 2);
    r = sqrt(1 - z^2);
    phi = mod(k,n) * increment;
    x = cos(phi) * r;
    y = sin(phi) * r;
    points(k,:) = [x,y,z];
end

if nargout > 1
    mindist = 1e6;
    for k1 = 1:n-1
        for k2 = k1+1:n
            dist = norm(points(k1,:)-points(k2,:));
            if dist < mindist
                mindist = dist;
            end
        end
    end
end