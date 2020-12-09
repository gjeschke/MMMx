function writeMRC(map,rez,filename,sizes,org)
% function WriteMRC(map,rez,filename)
% Write out a 2D image or a 3D volume as an MRC map file, for example for viewing in
% Chimera.  'map' is the 3D array, rez is the voxel size in angstroms.
%
% fixed for 2D output as well.  fs 28 Aug 07
% 
% % Test program: create an ellipsoid and store it as a map file.
%   [x y z]=ndgrid(-32:31);
%   m=(x.^2+(y*1.5).^2+z.^2)>20^2;
%   WriteMRC(m,1,'maptest.mrc');
%
% https://ch.mathworks.com/matlabcentral/fileexchange/50091-subspaceem-a-fast-maximum-a-posteriori-algorithm-for-cryo-em-single-particle-reconstruction

if nargin == 5
    f = WriteMRCHeader(map,rez,filename,sizes,org);
elseif nargin == 4
    f = WriteMRCHeader(map,rez,filename,sizes);
else
    f = WriteMRCHeader(map,rez,filename);
end
fwrite(f,map,'single');
fclose(f);

return;