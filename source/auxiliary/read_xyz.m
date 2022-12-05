function [ecoor,header] = read_xyz(fname,CB3D)
%
% Reads xyz file, also the Cartesian coordinate version as output by ChemBio3D
% output is empty if read fails
%
% fname   name of input file, if the file name has no extension, .xyz is
%         appended, for ChemBio3D files, .cc1 is appended
% CB3D    flag that indicates a ChemBio3D format (no header line), defaults
%         to false (standard xyz format)
%
% ecoor   extended coordinates, 1st column is the atom number, columns 2:4
%         are Cartesian coordinates in Angstroem
% header  header line as string
%
% (c) G. Jeschke, 2009-2013

ecoor = [];
header = '';

if ~exist('CB3D','var'), CB3D = false; end

if CB3D
    header = sprintf('Headerless ChemBio3D Cartesian coordinates %s',fname);
end

pse=' HHeLiBe B C N O FNeNaMgAlSi P SClAr KCaScTi VCrMnFeCoNiCuZnGaGeAsSeBrKrRbSr YZrNbMoTcRuRhPdAgCdInSnSbTe IXeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTa WReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPa UNpPuAmCmBkCfEsFmMdNoLr';

if ~contains(fname,'.')
    if CB3D
        fname=strcat(fname,'.cc1');
    else
        fname=strcat(fname,'.xyz');
    end
end
rfile = fopen(fname,'r');
if rfile == -1
    fprintf(2,'### ERROR: File could not be opened. ###\n');
    return;
end
    
line=fgets(rfile); % atom number line
if line==-1
    fprintf(2,'### ERROR: File shorter than expected. ###\n');
    fclose(rfile); return;
end
m=str2double(line);

if ~CB3D
    header = fgetl(rfile);
end

ecoor=zeros(m,4);
for kk=1:m
    line=fgets(rfile);
    if line==-1
        disp('### ERROR: Archive file shorter than expected. No output written. ###');
        fclose(rfile); return;
    end
    [symb,coordinates]=strtok(line);
    symb=deblank(symb);
    if length(symb)<2, symb=[' ' symb]; end
    element=(strfind(pse,symb)+1)/2;
    ecoor(kk,1)=element;
    aco=str2num(coordinates);
    if length(aco)~=3
        disp('### ERROR: File corrupted, wrong coordinate format. No output written. ###');
        fclose(rfile); return;
    end
    ecoor(kk,2:4)=aco;
end
