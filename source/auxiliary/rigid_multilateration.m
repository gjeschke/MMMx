function rigid_multilateration(restraints,logfid)
%
% fom = RIGID_MULTILATERATION(restraints,logfid)
%
% Locates a site by GPS-like multilateration, assuming a rigid polyhedron
% of reference points (beacons, like GPS satellites)
%
% restraints    restraints for localization, struct
%               .save_name  file name for output cube file 
%               .p_model    probability included by output isodensity
%                           surface
%               .grid       number of grid points, defaults to 176
%               .cube       cube size (Ansgtroem), defaults to 75 Angstroem
%               .reference  vector of k restraint blocks, each for N(k) 
%                           reference points, struct
%                           .r      vector (1,N(k)), mean distances between 
%                                   reference points and site (Angstroem)
%                           .sigr   vector (1,N(k)), standard deviations of 
%                                   reference-point-to-site distances
%                           .rax    cell (1,N(k)) of distance axes for full  
%                                   distributions
%                           .distr  cell (1,N(k)) of full distributions
%                           .coor   array(N(k),3) of Cartesian reference
%                                   point coordinates
% logfid        file identifier for log file
%
% G. Jeschke, 2021

% direct log to Matlab command window, if no file identifier was provided
if ~exist('logfid','var') || isempty(logfid)
    logfid = 1;
end

% set defaults
grid_points = 176;
cube_size = 75;

if isfield(restraints,'grid') && ~isempty(restraints.grid)
    grid_points = restraints.grid;
end
if isfield(restraints,'cube') && ~isempty(restraints.cube)
    cube_size = restraints.cube;
end

% reform restraints into single vectors/arrays
ref_points = zeros(1000,3);
distances = zeros(1000,1);
r_axes = cell(1000,1);
distributions = cell(1000,1);
np = 0;
for kb = 1:length(restraints.reference)
    for nr = 1:length(restraints.reference(kb).r)
        np = np + 1;
        ref_points(np,:) = restraints.reference(kb).coor{nr};
        distances(np) = restraints.reference(kb).r(nr);
        r_axes{np} = restraints.reference(kb).rax{nr};
        distributions{np} = restraints.reference(kb).distr{nr};
    end
end
ref_points = ref_points(1:np,:);
distances = distances(1:np);

% perform multilateration
[point, sos, singular] = multilaterate(ref_points,distances);

[m,~] = size(point);
if m == 2
    fprintf(logfid,'### Ambiguous location with only three reference points ###\n'); 
    fprintf(logfid,'Site was located either at mean coordinate (%5.1f,%5.1f,%5.1f) %s\n',point(1,:),char(197));
    fprintf(logfid,'                     or at mean coordinate (%5.1f,%5.1f,%5.1f) %s\n',point(2,:),char(197));
else
    fprintf(logfid,'Site was located at mean coordinate (%5.1f,%5.1f,%5.1f) %s\n',point,char(197));
    fprintf(logfid,'The inconsistency is %5.1f %s\n',sqrt(sos/np),char(197));
end
if singular
    fprintf(logfid,'### Warning: A singular matrix was encountered. Restraints are insufficient. ###\n');
end

% set grid size
nx = grid_points;
ny = nx;
nz = nx;

most_probable = [0,0,0];
for k=1:m
    max_probability = 0;
    %initialize cube
    cube = zeros(ny,nx,nz);
    x = linspace(point(k,1) - cube_size/2, point(k,1) + cube_size/2, nx);
    y = linspace(point(k,2) - cube_size/2, point(k,2) + cube_size/2, ny);
    z = linspace(point(k,3) - cube_size/2, point(k,3) + cube_size/2, nz);
    for kx = 1:nx
        for ky = 1:ny
            for kz = 1:nz
                coor = [x(kx),y(ky),z(kz)];
                diff = ref_points - reprowvector(coor,np);
                cdist = sqrt(sum(diff.^2,2));
                prob = get_probabilities(cdist,r_axes,distributions);
                cube(kx,ky,kz) = prod(prob);
                if cube(kx,ky,kz) > max_probability
                    max_probability = cube(kx,ky,kz);
                    most_probable = coor;
                end
            end
        end
    end
    cube = cube/max_probability; % normalize
    sdens=sum(sum(sum(cube)));
    % find the isodensity value that defines a fraction p_model of total
    % probability density
    level = restraints.p_model;
    level0 = level;
    for kl=1:99
        mask=(cube>=kl/100);
        test=sum(sum(sum(mask.*cube)));
        if test<=level0*sdens
            level = kl/100;
            break;
        end
    end
    [xg,yg,zg]=meshgrid(x,y,z);
    cnorm = sum(sum(sum(cube)));
    xm=sum(sum(sum(xg.*cube)))/cnorm;
    ym=sum(sum(sum(yg.*cube)))/cnorm;
    zm=sum(sum(sum(zg.*cube)))/cnorm;
    fprintf(logfid,'Center of gravity of density cloud %i is (%5.1f,%5.1f,%5.1f) %s\n',k,xm,ym,zm,char(197));
    fprintf(logfid,'A probability density isosurface at level %6.4f includes a fraction %4.2f of the total probability\n',level,restraints.p_model);
    fprintf(logfid,'The highest location probability is at (%4.1f, %4.1f, %4.1f) %s.\n',most_probable,char(197));
    m1 = max(max(cube(1,:,:)));
    m2 = max(max(cube(end,:,:)));
    m3 = max(max(cube(:,1,:)));
    m4 = max(max(cube(:,end,:)));
    m5 = max(max(cube(:,:,1)));
    m6 = max(max(cube(:,:,end)));
    mm = max(max(max(cube)));
    border_prop = max([m1,m2,m3,m4,m5,m6])/mm;
    fprintf(logfid,'The maximum relative probability at density cube border is %5.2f%%\n',100*border_prop);
    cube_name = restraints.save_name;
    [pname,fname,ext] = fileparts(cube_name);
    if m > 1        
        cube_name = fullfile(pname,sprintf('%s_%i%s',fname,k,ext));
    end
    switch ext
        case '.mrc'
            writeMRC(cube,x(2)-x(1),cube_name,[nx,ny,nz],[x(1),y(1),z(1)]);
        case '.mat'
            save(cube_name,'cube','x','y','z');
    end
end

function [point, sos, singular] = multilaterate(ref_points,dist)
% point = multilaterate(ref_points,distances)
%
% Determination of point coordinates by multilateration from three
% (trilateration) or more reference points and uncertain distances to these
% reference points
% 
% Input:
% ref_points    array [np,3] of Cartesian reference point coordinates
% dist          vector of the np distances to the reference points
% 
% Output:
% point         [1,3] row vector of Cartesian coordinates of the
%               multilaterated point, empty output for less than three
%               reference points, for three reference points, a [2,3] array
%               of two solutions or an empty array (no solution) may be
%               returned
% sos           sum of square distance errors in the multilateration case,
%               zero for successful trilateration, negative for unsuccesful
%               trilateration, zero for less than two reference points
% singular      flag that shows if a poorly conditioned matrix is
%               encountered in mutilateration (true) or not (false), this
%               happens for instance for multilateration from 4 points
%               within a plane, as far as I see this does not pose problems
%               for this particular task
%
% G. Jeschke, 25.4.2013

% the following uses the linearization of the problem by W. S. Murphy Jr.,
% Master thesis, Colorado School of Mines, 2007 for obtaining a starting
% coordinate for the standard non-linear least squares solution of the
% multilateration problem

testit = true;

sos = 0;
point = [];
singular = false;

[np,~] = size(ref_points);

if np<3
    return;
end

if np == 3 && ~testit % trilateration, actually based on Wikipedia article
    % go to standard frame
    origin = ref_points(1,:);
    coor = ref_points - repmat(origin,3,1);
    xe = coor(2,:)/norm(coor(2,:)); % unit vector along new x direction
    xye = coor(3,:)/norm(coor(3,:)); % another unit vector in xy plane
    if abs(sum(xe.*xye))>acos(5*pi/180) % points are almost collinear
        sos = -1;
        return;
    end
    ze = cross(xe,xye); % unit vector along new z direction
    ze = ze/norm(ze); % defensive programming
    ye = cross(ze,xe);
    ye = ye/norm(ye);
    Rp = [xe;ye;ze];
    coor1 = coor*Rp'; % coordinates in standard frame, point 1 origin, 
                      % point 2 on x axis, point 3 in xy plane
    d = coor1(2,1);
    i = coor1(3,1);
    j = coor1(3,2);
    x = (dist(1)^2 - dist(2)^2+d^2)/(2*d);
    y = (dist(1)^2 - dist(3)^2 + i^2 + j^2)/(2*j) - i*x/j;
    z = sqrt(dist(1)^2 - x^2 - y^2);
    if imag(z) > 10*eps % no solution
        sos = -2;
        return;
    elseif z < 10*eps % single solution
        point = [x,y,z]*Rp + origin;
    else % two solutions
        p1 = [x,y,real(z)]*Rp + origin;
        p2 = [x,y,-real(z)]*Rp + origin;
        point = [p1;p2];
    end
elseif np == 3 && testit
        % go to standard frame
    origin = ref_points(3,:);
    coor = ref_points - repmat(origin,3,1);
    xe = coor(1,:)/norm(coor(1,:)); % unit vector along new x direction
    xye = coor(2,:)/norm(coor(2,:)); % another unit vector in xy plane
    if abs(sum(xe.*xye))>acos(5*pi/180) % points are almost collinear
        sos = -1;
        return;
    end
    ze = cross(xe,xye); % unit vector along new z direction
    ze = ze/norm(ze); % normal vector to the plane
    A = zeros(np-1,3);
    b = zeros(np-1,1);
    for k=2:np
        A(k-1,:) = ref_points(k,:)-ref_points(1,:); % Eq. (2.9)
        dk1 = norm(ref_points(k,:)-ref_points(1,:)); % Eq. (2.4)
        b(k-1) = (dist(1)^2 - dist(k)^2 + dk1^2)/2; % Eq. (2.8)
    end

    if cond(A) < 1e10
        x = A\b;
    else
        x = pinv(A)*b;
        singular = true;
    end
    pt_0 = x' + ref_points(1,:);

    [point1, sos] = fminsearch(@sos_distances,pt_0,[],ref_points,dist);
    point2 = point1 + 2*dot(ref_points(1,:) - point1, ze)*ze;
    point = [point1;point2];
else % multilateration
    A = zeros(np-1,3);
    b = zeros(np-1,1);
    for k=2:np
        A(k-1,:) = ref_points(k,:)-ref_points(1,:); % Eq. (2.9)
        dk1 = norm(ref_points(k,:)-ref_points(1,:)); % Eq. (2.4)
        b(k-1) = (dist(1)^2 - dist(k)^2 + dk1^2)/2; % Eq. (2.8)
    end

    if cond(A) < 1e10
        x = A\b;
    else
        x = pinv(A)*b;
        singular = true;
    end
    pt_0 = x' + ref_points(1,:);

    [point, sos] = fminsearch(@sos_distances,pt_0,[],ref_points,dist);
end

function sos = sos_distances(pt,ref_points,distances)

[np,~] = size(ref_points);
diff = ref_points - repmat(pt,np,1);
dist_found = sqrt(sum(diff.^2,2));
sos = sum((dist_found-distances).^2);

function prob = get_probabilities(cdist,r_axes,distributions)

prob = zeros(1,length(cdist));
for k=1:length(cdist)
    [~,poi] = min(abs(r_axes{k} - cdist(k)));
    prob(k) = distributions{k}(poi);
end
