function fname = download_mmcif(ident)

fname = '';

query = sprintf('https://files.rcsb.org/download/%s.cif.gz',lower(ident));
fname0 = [lower(ident) '.cif.gz'];
try
    websave(fname0,query);
catch
    return;
end
try
    gunzip(fname0);
catch
    return;
end
fname = sprintf('%s.cif',lower(ident));
% delete the zipped file
try
    delete(fname0);
catch
    return
end