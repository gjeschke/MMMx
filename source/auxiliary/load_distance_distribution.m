function exp_data = load_distance_distribution(fname)
% loads distance distribution and converts distance axis to Angstroem, if
% it appears to be in nanometers

exp_data = load(fname);
r = exp_data(:,1);
% convert to Angstroem, if maximum distance is unusually short
if max(r) < 20
    exp_data(:,1) = 10*r;
end
