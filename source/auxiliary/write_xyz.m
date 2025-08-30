function write_xyz(fname,ecoor,comment)
% function write_xyz(fname,ecoor,comment)
%
% (c) G. Jeschke, 2008
%

pse=' HHeLiBe B C N O FNeNaMgAlSi P SClAr KCaScTi VCrMnFeCoNiCuZnGaGeAsSeBrKrRbSr YZrNbMoTcRuRhPdAgCdInSnSbTe IXeCsBaLaCePrNdPmSmEuGdTbDyHoErTmYbLuHfTa WReOsIrPtAuHgTlPbBiPoAtRnFrRaAcThPaUNpPuAmCmBkCfEsFmMdNoLr';


header=sprintf('Coordinates for %s',comment);

outname=strcat(fname,'.xyz');
[m,~]=size(ecoor);

wfile=fopen(outname,'w');
fprintf(wfile,'%i\n',m);
fprintf(wfile,'%s\n',header);
for k=1:m,
    z=ecoor(k,1);
    symb=pse(2*z-1:2*z);
    fprintf(wfile,'%3s%15.6f%14.6f%14.6f\n',symb,ecoor(k,2),ecoor(k,3),ecoor(k,4));
end;
fclose(wfile);