function [chi2,outname,status,result,fit] = fit_SAXS_by_crysol(datafile,pdbfile,options)
%

clear options
options.sm = 4.8016;
options.lm = 30;
options.fb = 18;

if ~exist('options','var') || ~isfield(options,'err')
    options.err = false;
end

poi = strfind(pdbfile,'.pdb');
if isempty(poi)
    outname = strcat(pdbfile,'00.fit');
    pdbfile = strcat(pdbfile,'.pdb');
else
    outname = [pdbfile(1:poi-1) '00.fit'];
end

s=which('crysol.exe');
cmd=[s ' ' pdbfile ' ' datafile ' -cst'];

if isfield(options,'sm')
    cmd = sprintf('%s -sm %6.3f',cmd,options.sm);
end

if isfield(options,'lm')
    cmd = sprintf('%s -lm %i',cmd,options.lm);
end

if isfield(options,'fb')
    cmd = sprintf('%s -fb %i',cmd,options.fb);
end

[status,result]=dos(cmd);

% chi = [];
% poi = strfind(result,'Chi:');
% if ~isempty(poi),
%     rem = textscan(result(poi+4:end),'%s');
%     args = rem{1};
%     chi = str2double(char(args(1)));
% end;
% 
% poi = strfind(result,'Data fit       saved to file');
% if ~isempty(poi),
%     rem = textscan(result(poi+length('Data fit       saved to file'):end),'%s');
%     args = rem{1};
%     outname = char(args(1));
% end;

fit = zeros(10000,4);

chi2 = [];

fid = fopen(outname);
if fid==-1
    fit = [];
    return;
end
nl=0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    %         fprintf(1,'%s\n',tline); % echo for debugging
    if nl > 0 % skip first line
        dataset = str2num(tline);
        ncol = length(dataset);
        fit(nl,1:ncol) = dataset;
    else
        poi = strfind(tline,'Chi^2:');
        if ~isempty(poi)
            rem = textscan(tline(poi+6:end),'%s');
            args = rem{1};
            chi2 = str2double(char(args(1)));
        else
            poi = strfind(tline,'Chi:');
            rem = textscan(tline(poi+4:end),'%s');
            args = rem{1};
            chi = str2double(char(args(1)));
            chi2 = chi^2;
        end
    end
    nl=nl+1;
end
fit = fit(1:nl-1,:);
poi = 1;
for k = 1:nl-1
    if fit(k,2) == 0 || fit(k,3) == 0
        poi = k;
    else
        break
    end
end
if ncol == 4 % remove error column from output of newer CRYSON
    if options.err
        fit = fit(:,[1 2 4 3]);
    else
        fit = fit(:,[1 2 4]);
    end
end
fit = fit(poi+1:end,:);
fclose(fid);

