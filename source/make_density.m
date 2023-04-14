function [x,y,z,cube] = make_density(entity,fname,address,resolution)
%
% MAKE_DENSITY Electron density map for an ensemble 
%
%   [x,y,z,cube] = MAKE_DENSITY(entity,fname)
%   Computes density and stores it in MRC cube file fname
%
%   [x,y,z,cube] = MAKE_DENSITY(entity,fname,address)
%   Restricts computation to the part selected by address
%
% INPUT
% entity       entity in an MMMx format, must be provided
% fname        file name for output, extension .mrc is
%              appended, if none is present
% address      MMMx address, defaults to 'everything'
% resolution   (optional) grid resolution in Angstroem
%
% OUTPUT
% x,y,z        cube axes (Angstrom)
% cube         density cube
%
% algorithm according to:
% R. Briones, C. Blau, C. Kutzner, B. L. de Groot, C. Aponte-Santamarı,
% Biophys. J. 116, 4–11 (2019). https://doi.org/10.1016/j.bpj.2018.11.3126
%
% Gaussian coefficients for point-spread functions are taken from the SI of
% the same paper and are provided via function point_spread.m
% all atoms except for H, N, O, P, and S are treated as C

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke


% default selection is everything (all chains)
if ~exist('address','var') || isempty(address)
    address = '(*)';
end

% append extension .mrc to file name, if there is none
if ~contains(fname,'.')
    fname = [fname '.mrc'];
end

slater_generic = 0.75; % generic Slater radius for heavy atoms, 75 pm to include hydrogen

rez = slater_generic/2; % default resolution, very high
if exist('resolution','var') && ~isempty(resolution)
    rez = resolution;
end

C = length(entity.populations); % number of conformers
min_xyz = [1e6,1e6,1e6];
max_xyz = [-1e6,-1e6,-1e6];
for c = 1:C
    coor = get_coor(entity,sprintf('{%i}%s',c,address),true);
    max_xyz0 = max(coor) + 3*slater_generic;
    min_xyz0 = min(coor) - 3*slater_generic;
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
cube = zeros(nx,ny,nz);
x = linspace(min_xyz(1), max_xyz(1), nx);
y = linspace(min_xyz(2) , max_xyz(2), ny);
z = linspace(min_xyz(3), max_xyz(3), nz);

for c = 1:C
    fprintf(1,'Conformer %i\n',c);
    pop = entity.populations(c);
    [coor,~,elements] = get_coor(entity,sprintf('{%i}%s',c,address),true);
    [n,~] = size(coor);
    for at = 1:n
        AB = point_spread(elements(at));
        xpsf = zeros(size(x));
        ypsf = zeros(size(y));
        zpsf = zeros(size(z));
        for g = 1:4
            xpsf = xpsf + AB(g,1)*exp(-AB(g,2)*(x-coor(at,1)).^2); % point-spread function along x
            ypsf = ypsf + AB(g,1)*exp(-AB(g,2)*(y-coor(at,2)).^2); % point-spread function along y
            zpsf = zpsf + AB(g,1)*exp(-AB(g,2)*(z-coor(at,3)).^2); % point-spread function along z
        end
        transverse = xpsf.'*ypsf; % transverse point-spread function
        for kz = 1:nz
            cube(:,:,kz) = cube(:,:,kz) + pop*zpsf(kz)*transverse;
        end
    end
end

[~,~,ext] = fileparts(fname);
switch ext
    case '.mrc'
        writeMRC(cube,rez,fname,[nx,ny,nz],[x(1),y(1),z(1)]);
    case '.mat'
        save(fname,'cube','x','y','z');
end