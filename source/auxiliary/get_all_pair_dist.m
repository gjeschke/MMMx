function pair_dist = get_all_pair_dist(a,b)
% pair_dist = get_all_pair_dist(a,b)
%
% get all distances between atom pairs for two coordinate arrays a and b
%
% code provided by Stefan Stoll

[m1,~] = size(a); % get sizes of the coordinates arrays
[m2,~] = size(b);

a2 = repmat(sum(a.^2,2),1,m2);
b2 = repmat(sum(b.^2,2),1,m1).';
pair_dist = sqrt(abs(a2 + b2 - 2*a*b.'));
