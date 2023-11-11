function [disorder,CheZOD] = get_disorder_SETH(database,UniProtID,threshold)

if ~exist('threshold','var') || isempty(threshold)
    threshold = 8;
end

CheZOD = [];
fid = fopen(database,'rt');
while 1
    tline = fgetl(fid);
    if ~ischar(tline)
        break
    end
    args = split(tline,'|');
    if length(args) > 1
        if strcmpi(args{2},UniProtID)
            tline = fgetl(fid);
            CheZOD = str2num(tline); %#ok<ST2NM>
        end
    end
end
if ~isempty(CheZOD)
    disorder = CheZOD < threshold;
end
fclose(fid);