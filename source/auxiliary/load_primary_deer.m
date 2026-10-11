function exp_data = load_primary_deer(fname)
% loads primary DEER data

[~,~,ext] = fileparts(fname);

switch ext
    case '.csv'
        exp_data = read_my_csv(fname);
    case '.dat'
        exp_data = load(fname);
end


function exp_data = read_my_csv(fname)
data = zeros(10000,5);
dpoi = 0;
format = 'coventional';
fid = fopen(fname);
while 1
    tline = fgetl(fid);
    if ~ischar(tline) 
        break 
    end
    if contains(tline,'%') % skip header line of DeerAnalysis2022 output
        continue
    end
    args = split(tline,',');
    if isnan(str2double(args{1})) % recognize DeerAnalysis2026 output and skip header line
        format = 'DeerAnalysis2026';
        continue
    end
    dpoi = dpoi + 1;
    if strcmpi(format,'conventional')
        for k = 1:4
            data(dpoi,k) = str2double(args{k});
        end
    elseif strcmpi(format,'DeerAnalysis2026')
        for k = 1:5
            data(dpoi,k) = str2double(args{k});
        end
    end
end
fclose(fid);
if strcmpi(format,'DeerAnalysis2026')
   [~,idx0] = max(data(:,4));
   data(:,2) = data(:,2)/data(idx0,2);
   data(:,4) = data(:,5);
   data(:,1) = data(:,1) - data(idx0,1);
   exp_data = data(idx0:dpoi,1:4);
else
   exp_data = data(1:dpoi,1:4);
end
