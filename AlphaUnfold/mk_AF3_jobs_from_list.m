function mk_AF3_jobs_from_list(shortlist,domain)

number = 50000;
UPID = cell(1,number);

fid = fopen(shortlist);
fgetl(fid);

matches = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    matches = matches + 1;
    args = split(tline,',');
    UPID{matches} = args{1};
end

fclose(fid);

UPID = UPID(1:matches);

if ~exist(domain,'dir')
    mkdir(domain);
end
cd(domain);

done = 0;
for n = 1:matches
    mk_job(UPID{n},n);
    done = done + 1;
    if mod(done,100) == 0
        fprintf(1,'%i of %i AF3 jobs prepared\n',done,matches);
    end
end

cd('..');

fprintf(1,'%i AF3 jobs were prepared\n',matches);

function mk_job(UniProtID,n)

seq = get_sequence(UniProtID);
seed = round(1e9*rand);
fname = sprintf('J%i_%s.json',n,UniProtID);
fid = fopen(fname,'wt');
fprintf(fid,'[{\n');
fprintf(fid,' "name": "J%i_%s",\n',n,UniProtID);
fprintf(fid,' "modelSeeds": [\n');
fprintf(fid,'  "%i"\n',seed);
fprintf(fid,' ],\n');
fprintf(fid,' "sequences": [\n');
fprintf(fid,'   {\n');
fprintf(fid,'    "proteinChain": {\n');
fprintf(fid,'    "sequence":\n');
fprintf(fid,'    "%s",\n',seq);
fprintf(fid,'    "count": 1\n');
fprintf(fid,'   }\n');
fprintf(fid,'  }');
fprintf(fid,'\n');
fprintf(fid,' ],\n');
fprintf(fid,'"dialect": "alphafoldserver",\n');
fprintf(fid,'"version": 1\n');
fprintf(fid,'}\n');
fprintf(fid,']\n');
fclose(fid);
