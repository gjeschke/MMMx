function properties = get_orca_properties(OrcaFileName)
% get_orca_properties   Import a subset of ORCA properties into a struct
% variable
%
%  properties = get_orca_properties(OrcaFileName)
%  Sys = orca2easyspin(OrcaFileName,HyperfineCutoff)
%
%  Loads the most important parameters from an ORCA property text
%  output file given in OrcaFileName and returns them as a Matlab
%  structure
%

properties.energy = 0;

token{1} = 'SCF_Energy';
token{2} = 'EPRNMR_GTensor';
token{3} = 'EPRNMR_ATensor';
token{4} = 'EPRNMR_QTensor';

rfid = fopen(OrcaFileName);
if rfid == -1
    properties = [];
    return;
end
while 1
    tline = fgetl(rfid);
    if ~ischar(tline) 
        break 
    end
    if tline(1) == '$' % this is a token line
        ctoken = strip(tline(2:end));
        switch ctoken
            case token{1}
                properties = get_SCF_energy(properties,rfid);
            case token{2}
                properties = get_g_tensor(properties,rfid);
            case token{3}
                properties = get_A_tensors(properties,rfid);
            case token{4}
                properties = get_Q_tensors(properties,rfid);
            otherwise
                fprintf(1,'Token %s is skipped\n',ctoken);
        end
    end
end
fclose(rfid);

function properties = get_SCF_energy(properties,rfid)

key = 'SCF Energy:';
key_found = false;
while ~key_found
    tline = fgetl(rfid);
    if ~ischar(tline) 
        break 
    end
    poi = strfind(tline,key);
    if ~isempty(poi)
        key_found = true;
        properties.energy = str2double(tline(poi+length(key):end));
    end
end

function properties = get_g_tensor(properties,rfid)

key = 'Spin multiplicity:';
key_found = false;
while ~key_found
    tline = fgetl(rfid);
    if ~ischar(tline) 
        break 
    end
    poi = strfind(tline,key);
    if ~isempty(poi)
        key_found = true;
        properties.spin = (str2double(tline(poi+length(key):end))-1)/2;
    end
end
fgetl(rfid);
fgetl(rfid);
properties.g_tensor = zeros(3);
for k1 = 1:3
    tline = fgetl(rfid);
    args = split(tline);
    for k2 = 1:3
        properties.g_tensor(k1,k2) = str2double(args{2+k2});
    end
end
g_symm = (properties.g_tensor+properties.g_tensor')/2;
[evecs,eval] = eig(g_symm);
properties.g_values = diag(eval);
properties.g_frame = evecs;

function properties = get_A_tensors(properties,rfid)

key = 'Number of stored nuclei';
key_found = false;
while ~key_found
    tline = fgetl(rfid);
    if ~ischar(tline) 
        break 
    end
    poi = strfind(tline,key);
    if ~isempty(poi)
        key_found = true;
        nuclei = round(str2double(tline(poi+length(key):end)));
    end
end
properties.A_tensors = zeros(3*nuclei,3);
properties.A_frames = zeros(3*nuclei,3);
properties.A_values = zeros(nuclei,3);
properties.nuclei = cell(nuclei,1);
properties.nuc_spin = zeros(nuclei,1);
properties.isotopes = zeros(nuclei,1);
properties.input = zeros(nuclei,1);
fgetl(rfid);
fgetl(rfid);
for n = 1:nuclei
    pointer = 3*(n-1);
    tline = fgetl(rfid);
    args = split(tline);
    properties.input(n) = str2double(args{3});
    properties.nuclei{n} = args{4};
    tline = fgetl(rfid);
    args = split(tline);
    properties.isotopes(n) = str2double(args{3});
    tline = fgetl(rfid);
    args = split(tline);
    properties.nuc_spin(n) = str2double(args{4});
    for k = 1:3
        fgetl(rfid);
    end
    A_tensor = zeros(3);
    for k1 = 1:3
        tline = fgetl(rfid);
        args = split(tline);
        for k2 = 1:3
            A_tensor(k1,k2) = str2double(args{2+k2});
        end
    end
    properties.A_tensors(pointer+1:pointer+3,:) = A_tensor;
    A_symm = (A_tensor+A_tensor')/2;
    [evecs,eval] = eig(A_symm);
    properties.A_values(n,:) = diag(eval);
    properties.A_frames(pointer+1:pointer+3,:) = evecs;
    for k = 1:9 % skip remaining lines with redundant information
        fgetl(rfid);
    end
end
        
function properties = get_Q_tensors(properties,rfid)

key = 'Number of stored nuclei';
key_found = false;
while ~key_found
    tline = fgetl(rfid);
    if ~ischar(tline) 
        break 
    end
    poi = strfind(tline,key);
    if ~isempty(poi)
        key_found = true;
        nuclei = round(str2double(tline(poi+length(key):end)));
    end
end
properties.Q_tensors = zeros(3*nuclei,3);
properties.Q_frames = zeros(3*nuclei,3);
properties.Q_values = zeros(nuclei,3);
properties.Q_nuclei = cell(nuclei,1);
properties.Q_nuc_spin = zeros(nuclei,1);
properties.Q_isotopes = zeros(nuclei,1);
properties.Q_input = zeros(nuclei,1);
fgetl(rfid);
fgetl(rfid);
fgetl(rfid);
for n = 1:nuclei
    pointer = 3*(n-1);
    tline = fgetl(rfid);
    args = split(tline);
    properties.Q_input(n) = str2double(args{3});
    properties.Q_nuclei{n} = args{4};
    tline = fgetl(rfid);
    args = split(tline);
    properties.Q_isotopes(n) = str2double(args{3});
    tline = fgetl(rfid);
    args = split(tline);
    properties.Q_nuc_spin(n) = str2double(args{4});
    for k = 1:3
        fgetl(rfid);
    end
    Q_tensor = zeros(3);
    for k1 = 1:3
        tline = fgetl(rfid);
        args = split(tline);
        for k2 = 1:3
            Q_tensor(k1,k2) = str2double(args{2+k2});
        end
    end
    properties.Q_tensors(pointer+1:pointer+3,:) = Q_tensor;
    Q_symm = (Q_tensor+Q_tensor')/2;
    [evecs,eval] = eig(Q_symm);
    properties.Q_values(n,:) = diag(eval);
    properties.Q_frames(pointer+1:pointer+3,:) = evecs;
    for k = 1:9 % skip remaining lines with redundant information
        fgetl(rfid);
    end
end
        

