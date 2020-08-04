function dssp = get_dssp(fname)

dssp=[];

dospath = which('dssp-2.0.4-win32.exe');
if isempty(dospath)
    dospath = which('dssp.exe');
end
if isempty(dospath)
    dospath = which('dsspcmbi.exe');
end
if ~isempty(dospath) % suppress this if DSSP not known
    infile = which(fname);
    poi = strfind(infile,'.pdb');
    if isempty(poi)
        poi = strfind(infile,'.ent');
    end
    if isempty(poi) || poi(1)<2
        basname = infile;
    else
        basname = infile(1:poi(1)-1);
    end
    dssp_file = [basname '.dssp'];
    cmd=[dospath ' ' infile ' ' dssp_file];
    s = dos(cmd);
    if s~=0 % dssp was not successful
        return
    else
        dssp = rd_dssp(dssp_file);
    end
end 


