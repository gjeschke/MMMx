function residues = rd_dssp(fname)
% function residues = rd_dssp(fname)
%
% Reads a DSSP (Kabsch/Sanders) file and stores the results in a structure
% for the residues
%
% fname     name of the DSSP file
% residues  results, an array of structures with the following fields:
%           .tag    residue tag (number and possibly insertion code
%           .chain  chain identifier
%           .slc    single-letter code
%           .sec    secondary structure assignment, a single character
%                   space (blank)   loop or irregular
%                   H               alpha helix
%                   B               isolated beta bridge
%                   E               extended strand
%                   G               3/10-helix
%                   I               pi helix
%                   T               hydrogen-bonded turn
%                   S               bend
%           .bp     two integer numbers of bridge partners
%           .sheet  sheet identifier
%           .acc    water-exposed surface in A^2
%           .NHO    N-H-->O hydrogen bonds (two) with subfields
%                   .hp     hydrogen bonding partner
%                   .energy energies (kcal/mol)
%           .OHN    O-->H-N hydrogen bonds (two) with subfields
%                   .hp     hydrogen bonding partner
%                   .energy energies (kcal/mol)
%           .tco    cosine of angle between C=O of residue i and residue
%                   i-1, about 1 for alpha helix, about -1 for beta-sheet
%           .kappa  virtual bond angle between C_alpha of residues
%                   i-2,i,i+2 (for bend definition)
%           .alpha  virtual torsion angle between four C_alpha atoms of
%                   residues i-1,i,i+1,i+2
%           .phi    backbone dihedral phi
%           .psi    backbone dihedral psi
%           
% see: http://swift.cmbi.kun.nl/gv/dssp/ for DSSP license and download
% citation: 
%           W. Kabsch, C. Sanders, Biopolymers, 1983, 22(12):2577-637.
%           Dictionary of protein secondary structure: pattern recognition 
%           of hydrogen-bonded and geometrical features.
%
% G. Jeschke, 2009

residues = [];

fid = fopen(fname,'r');

poi=0;
residues(5000).tag = '   ';
if fid~=-1
    isheader=1;
    while isheader
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        if contains(tline,'  #  RESIDUE')
            isheader=0;
        end
    end
    while ~isheader% residue information is only read if keyword was found
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        if length(tline)>=115
            poi=poi+1;
            if strcmp(tline(12),' '), tline(12) = 'Z'; end % patch missing chain ID
            residues(poi).tag = strtrim(tline(6:11));
            residues(poi).chain = tline(12);
            residues(poi).slc = tline(14);
            residues(poi).sec = tline(17);
            residues(poi).bp(1) = sscanf(tline(26:29),'%i');
            residues(poi).bp(2) = sscanf(tline(30:33),'%i');
            residues(poi).sheet = tline(34);
            residues(poi).acc = sscanf(tline(35:38),'%i');
            residues(poi).NHO(1).hp = sscanf(tline(40:45),'%i');
            residues(poi).NHO(1).energy = sscanf(tline(47:50),'%f');
            residues(poi).OHN(1).hp = sscanf(tline(51:56),'%i');
            residues(poi).OHN(1).energy = sscanf(tline(58:61),'%f');
            residues(poi).NHO(2).hp = sscanf(tline(62:67),'%i');
            residues(poi).NHO(2).energy = sscanf(tline(69:72),'%f');
            residues(poi).OHN(2).hp = sscanf(tline(73:78),'%i');
            residues(poi).OHN(2).energy = sscanf(tline(80:83),'%f');
            residues(poi).tco = sscanf(tline(84:91),'%f');
            residues(poi).kappa = sscanf(tline(92:97),'%f');
            residues(poi).alpha = sscanf(tline(98:103),'%f');
            residues(poi).phi = sscanf(tline(104:109),'%f');
            residues(poi).psi = sscanf(tline(110:115),'%f');
        end
    end
end
residues = residues(1:poi);
fclose(fid);