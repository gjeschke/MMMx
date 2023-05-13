function coor = rd_accutar(fname)
% RD_ACCUTAR reads output of the hydration prediction server Accutar
%
% coor = rd_accutar(fname)
%   provides coordinates of O atoms of generated water molecules
%
% Input:
%
% fname     file name of the PDB output file of Accutar, use
%           https://oa.accutarbio.com to generate it, extension .pdb is
%           appended if missing
%
% Output:
%
% coor      (n,3) coordinate array of n water oxygen atom coordinates (Å)
% 

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

[~,~,ext] = fileparts(fname);

% load file from PDB server if required
if isempty(ext)
    fname = strcat(fname,'.pdb');
end

fid = fopen(fname);

coor = zeros(10000,3);
nO = 0;

while 1
    tline = fgetl(fid);
    if ~ischar(tline) 
        break 
    end
    if length(tline) >= 54 && strcmpi(tline(1:6),'HETATM')
        if strcmpi(tline(18:20),'HOH') && tline(14) == 'O'
            nO = nO + 1;
            coor(nO,1) = str2double(tline(31:38));
            coor(nO,2) = str2double(tline(39:46));
            coor(nO,3) = str2double(tline(47:54));
        end
    end
end

coor = coor(1:nO,:);
