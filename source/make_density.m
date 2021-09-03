function make_density(entity,fname,address,resolution)
%
% MAKE_DENSITY Pseudo-electron density for an ensemble 
%
%   MAKE_DENSITY(entity,fname)
%   Computes density and stores it in MRC cube file fname
%
%   MAKE_DENSITY(entity,fname,address)
%   Restricts computation to the part selected by address
%
% INPUT
% entity       entity in an MMMx format, must be provided
% fname        file name for output, extension .mrc is
%              appended, if none is present
% address      MMMx address, defaults to 'everything'
% resolution   (optional) grid resolutiin in Angstroem
%

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

sig = 3/sqrt(8*log(2)); % 3 Angstr?m FWHH resolution

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
    coor = get_coor(entity,sprintf('{%i}%s',c,address),true);
    [n,~] = size(coor);
    for at = 1:n
        xpsf = exp(-(x-coor(at,1)).^2/(2*sig^2)); % point-spread function along x
        ypsf = exp(-(y-coor(at,2)).^2/(2*sig^2)); % point-spread function along y
        zpsf = exp(-(z-coor(at,3)).^2/(2*sig^2)); % point-spread function along y
        transverse = xpsf.'*ypsf; % transverse point-spread function
        for kz = 1:nz
            cube(:,:,kz) = cube(:,:,kz) + pop*zpsf(kz)*transverse;
        end
    end
end

writeMRC(cube,rez,fname,[nx,ny,nz],[x(1),y(1),z(1)]);