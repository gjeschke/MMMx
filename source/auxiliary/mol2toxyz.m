function mol2toxyz(fname)
% MOL2TOXYZ
%
%  mol2toxyz(fname)
%
% Converts a Sybyl .mol2 file to a .xyz coordinate file

[pname,basname,ext] = fileparts(fname);
if isempty(ext)
    fname = strcat(fname,'.mol2');
end
oname = fullfile(pname,sprintf('%s.xyz',basname));

ifid = fopen(fname,'rt');
ofid = fopen(oname,'wt');
fprintf(ofid,'Mol2 File %s converted to xyz\n',basname);

read_atoms = false;
while ~read_atoms
    tline = fgetl(ifid);
    if ~ischar(tline) 
        break 
    end
    if length(tline) >= 17 && strcmpi(tline(1:17),'@<TRIPOS>MOLECULE')
        fgetl(ifid);
        tline = fgetl(ifid);
        args = split(tline);
        atoms = str2double(args{1});
        fprintf(ofid,'%i\n',atoms);
    end
    if length(tline) >= 13 && strcmpi(tline(1:13),'@<TRIPOS>ATOM')
        read_atoms = true;
        for k = 1:atoms
            tline = fgetl(ifid);
            args = split(tline);
            while isempty(args{1})
                args = args(2:end);
            end
            x = str2double(args{3});
            y = str2double(args{4});
            z = str2double(args{5});
            fprintf(ofid,'%3s%15.6f%14.6f%14.6f\n',args{2},x,y,z);
        end
    end
end
fclose(ifid);
fclose(ofid);
