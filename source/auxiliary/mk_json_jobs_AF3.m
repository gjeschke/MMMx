function mk_json_jobs_AF3(fbasname,sequence,type,copies,seeds)

nseeds = length(seeds);
seeds = zeros(nseeds,1);
seeds(1) = round(1e9*rand);
for s = 2:nseeds
    seed = 0;
    unique = min(abs(seeds(1:s-1)-seed));
    while ~unique || seed < 1 || seed > 999999999
        seed = round(1e9*rand);
    end
    seeds(s) = seed;
end

for s = 1:nseeds
    fname = sprintf('%s_%i.json',fbasname,s);
    fid = fopen(fname,'wt');
    fprintf(fid,'[{\n');
    fprintf(fid,' "name": "%s_%i",\n',fbasname,s);
    fprintf(fid,' "modelSeeds": [\n');
    fprintf(fid,'  "%i"\n',seeds(s));
    fprintf(fid,' ],\n');
    fprintf(fid,' "sequences": [\n');
    for seq = 1:length(sequence)
        fprintf(fid,'   {\n');
        fprintf(fid,'    "%s": {\n',type{seq});
        fprintf(fid,'    "sequence":\n');
        fprintf(fid,'    "%s",\n',sequence{seq});
        fprintf(fid,'    "count": %i\n',copies{seq});
        fprintf(fid,'   }\n');
        fprintf(fid,'  }');
        if s < length(sequence)
            fprintf(fid,',');
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,' ]\n');
    fprintf(fid,'}\n');
    fprintf(fid,']\n');
    fclose(fid);
end