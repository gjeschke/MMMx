function [disorder,probability] = get_disorder_espritz(UniProtID,threshold)

fname = sprintf('sp%s.fasta.espritz',UniProtID);
fid = fopen(fname,'rt');
if fid < 0
    fname = sprintf('tr%s.fasta.espritz',UniProtID);
    fid = fopen(fname,'rt');
end
if fid < 0
    disorder = [];
    probability = [];
    return
end
disorder = zeros(10000,1);
probability = zeros(10000,1);
poi = 0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    poi = poi + 1;
    args = split(tline);
    if strcmpi(args{1},'D')
        disorder(poi) = 1;
    end
    probability(poi) = str2double(args{2});
end
disorder = disorder(1:poi);
probability = probability(1:poi);
if exist('threshold','var') && ~isempty(threshold)
    disorder = probability > threshold;
end
fclose(fid);