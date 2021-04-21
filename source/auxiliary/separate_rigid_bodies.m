function entity = separate_rigid_bodies(entity,rigid_bodies)
% function SEPARATE_RIGID_BODIES(entity,rigid_bodies)
%
% Separates the rigid bodies in an entity as preparation for a Rigi run
% this is needed in order to avoid erroneous clashes in in silico spin
% labelling of individual rigid bodies
%
% Input:
%
% entity     MMMx entity
% 
%
% output file name is determined by the # PDB command in the restraint
% file
%
% G. Jeschke, 10.9.2019


maxr = 0;
rb_poi = length(rigid_bodies);
centers = zeros(rb_poi,3);  % center coordinates of each rigid body

% determine the center coordinate of each rigid body
for krb = 1:rb_poi
    atnum = 0;
    bcenter = zeros(1,3);
    for kc = 1:length(rigid_bodies(krb).chains)
        address = sprintf('(%s)',rigid_bodies(krb).chains(kc));
        xyz = get_coor(entity,address);
        [m,~] = size(xyz);
        atnum = atnum + m;
        bcenter = bcenter + m*mean(xyz);
    end
    centers(krb,:) = bcenter/atnum;
end

% determine the extension of each rigid body with respect to its center
% coordinate
for krb = 1:rb_poi
    for kc = 1:length(rigid_bodies(krb).chains)
        address = sprintf('(%s)',rigid_bodies(krb).chains(kc));
        xyz = get_coor(entity,address);
        [m,~] = size(xyz);
        xyzs = xyz - repmat(centers(krb,:),m,1);
        rmax = max(sqrt(sum(xyzs.^2,2)));
        if rmax > maxr
            maxr = rmax;
        end
    end
end

% distribute rigid body midpoints with a fibonacci sphere 
[midpoints,mindist] = fibonacci_sphere(rb_poi);
midpoints = midpoints*2*maxr/mindist;

for krb = 1:rb_poi
    % determine shift of this rigid body
    shift = midpoints(krb,:) - centers(krb,:);
    for kc = 1:length(rigid_bodies(krb).chains)
        address = sprintf('(%s)',rigid_bodies(krb).chains(kc));
        [xyz,indices] = get_coor(entity,address);
        [m,~] = size(xyz);
        % shift the coordinates
        xyz = xyz + repmat(shift,m,1);
        entity.xyz(indices,:) = xyz;
    end
end
