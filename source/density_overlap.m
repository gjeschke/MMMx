function overlap = density_overlap(entity1,entity2,address,resolution)
%
% DENSITY_OVERLAP Computes pseudo-electron density overlap of two ensembles 
%
%   overlap = DENSITY_OVERLAP(entity1,entity2)
%   Computes density overlap of whole structures
%
%   overlap = DENSITY_OVERLAP(entity1,entity2,address)
%   Restricts computation to the part selected by address
%
%   overlap = DENSITY_OVERLAP(entity1,entity2,address,resolution)
%   Specifies resolution, default ist 1.5 Angstroem
%
% INPUT
% entity       entity in an MMMx format, must be provided
% fname        file name for output, extension .mrc is
%              appended, if none is present
% address      MMMx address, defaults to 'everything'
% resolution   (optional) grid resolution in Angstroem, defaults to 1.274
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke


% default selection is everything (all chains)
if ~exist('address','var') || isempty(address)
    address = '(*)';
end

sig = 3/sqrt(8*log(2)); % 3 Angstr?m FWHH resolution
rez = sig;

if exist('resolution','var') && ~isempty(resolution)
    rez = resolution;
end

C1 = length(entity1.populations); % number of conformers
C2 = length(entity2.populations); % number of conformers
min_xyz = [1e6,1e6,1e6];
max_xyz = [-1e6,-1e6,-1e6];
for c = 1:C1
    coor = get_coor(entity1,sprintf('{%i}%s',c,address),true);
    max_xyz0 = max(coor) + 3*sig;
    min_xyz0 = min(coor) - 3*sig;
    for k = 1:3
        if max_xyz(k) < max_xyz0(k)
            max_xyz(k) = max_xyz0(k);
        end
        if min_xyz(k) > min_xyz0(k)
            min_xyz(k) = min_xyz0(k);
        end
    end
end
for c = 1:C2
    coor = get_coor(entity2,sprintf('{%i}%s',c,address),true);
    max_xyz0 = max(coor) + 3*sig;
    min_xyz0 = min(coor) - 3*sig;
    for k = 1:3
        if max_xyz(k) < max_xyz0(k)
            max_xyz(k) = max_xyz0(k);
        end
        if min_xyz(k) > min_xyz0(k)
            min_xyz(k) = min_xyz0(k);
        end
    end
end
nx = round((max_xyz(1) - min_xyz(1))/rez);
ny = round((max_xyz(2) - min_xyz(2))/rez);
nz = round((max_xyz(3) - min_xyz(3))/rez);
cube1 = zeros(nx,ny,nz);
cube2 = zeros(nx,ny,nz);
x = linspace(min_xyz(1), max_xyz(1), nx);
y = linspace(min_xyz(2) , max_xyz(2), ny);
z = linspace(min_xyz(3), max_xyz(3), nz);

for c = 1:C1
    pop = entity1.populations(c);
    coor = get_coor(entity1,sprintf('{%i}%s',c,address),true);
    [n,~] = size(coor);
    for at = 1:n
        xpsf = exp(-(x-coor(at,1)).^2/(2*sig^2)); % point-spread function along x
        ypsf = exp(-(y-coor(at,2)).^2/(2*sig^2)); % point-spread function along y
        zpsf = exp(-(z-coor(at,3)).^2/(2*sig^2)); % point-spread function along y
        transverse = xpsf.'*ypsf; % transverse point-spread function
        for kz = 1:nz
            cube1(:,:,kz) = cube1(:,:,kz) + pop*zpsf(kz)*transverse;
        end
    end
end
cube1 = cube1/sum(sum(sum(cube1)));

for c = 1:C2
    pop = entity2.populations(c);
    coor = get_coor(entity2,sprintf('{%i}%s',c,address),true);
    [n,~] = size(coor);
    for at = 1:n
        xpsf = exp(-(x-coor(at,1)).^2/(2*sig^2)); % point-spread function along x
        ypsf = exp(-(y-coor(at,2)).^2/(2*sig^2)); % point-spread function along y
        zpsf = exp(-(z-coor(at,3)).^2/(2*sig^2)); % point-spread function along y
        transverse = xpsf.'*ypsf; % transverse point-spread function
        for kz = 1:nz
            cube2(:,:,kz) = cube2(:,:,kz) + pop*zpsf(kz)*transverse;
        end
    end
end
cube2 = cube2/sum(sum(sum(cube2)));

overlap = sum(sum(sum(min(cube1,cube2))));
