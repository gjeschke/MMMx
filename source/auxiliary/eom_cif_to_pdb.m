function entity = eom_cif_to_pdb(fname)

[filepath,filename,ext] = fileparts(fname);
if strcmpi(ext,'.cif')
    iname = fname;
else
    iname = strcat(fullfile(filepath,filename),'.cif');
end
oname = strcat(fullfile(filepath,filename),'.pdb');

fid = fopen(iname,'rt');
ofid = fopen(oname,'wt');

while 1
    if fid~=-1
        tline = fgetl(fid);
        if ~ischar(tline), break, end
        if length(tline) >= 6
            if strcmpi(tline(1:4),'ATOM') % || strcmpi(tline(1:6),'HETATM') to be implemented later
                args = split(tline);
                anum = str2double(args{2});
                elm = args{3};
                atag = args{4};
                rtag = args{6};
                ctag = args{7};
                rnum = str2double(args{8});
                x = str2double(args{10});
                y = str2double(args{11});
                z = str2double(args{12});
                occupancy = str2double(args{13});
                Bfactor = str2double(args{14});
                nline = sprintf('ATOM  %5i  %s',anum,atag);
                nline = pad(nline,17,'right');
                nline = sprintf('%s%s %s%4i%12.3f%8.3f%8.3f%6.2f%6.2f',nline,rtag,ctag,rnum,x,y,z,occupancy,Bfactor);
                nline = pad(nline,76,'right');
                nline = sprintf('%2s',nline,elm);
                fprintf(ofid,'%s\n',nline);
            end
        end
    else
        break
    end
end

fclose(fid);
fclose(ofid);

if nargout > 0
    entity = get_pdb(oname);
end
