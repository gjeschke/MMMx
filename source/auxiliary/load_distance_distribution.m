function exp_data = load_distance_distribution(fname)
% loads distance distribution and converts distance axis to Angstroem, if
% it appears to be in nanometers

exp_data = read_my_csv(fname);
r = exp_data(:,1);
% convert to Angstroem, if maximum distance is unusually short
if max(r) < 20
    exp_data(:,1) = 10*r;
end

function data = read_my_csv(fname)
data = zeros(10000,4);
dpoi = 0;
fid = fopen(fname);
while 1
    tline = fgetl(fid);
    if ~ischar(tline) 
        break 
    end
    poi = strfind(tline,'%');
    if ~isempty(poi)
        continue
    end
    args = split(tline,',');
    dpoi = dpoi + 1;
    for k = 1:4
        data(dpoi,k) = str2double(args{k});
    end
end
fclose(fid);
data = data(1:dpoi,:);