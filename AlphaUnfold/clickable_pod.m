[filename, pathname] = uigetfile('_disorder_characteristics.csv', 'Select file');
if isequal(filename,0)
    disp('User selected Cancel');
else
    fullpath = fullfile(pathname, filename);
end

data = get_csv(fullpath,'data.cell');
[proteins,~] = size(data);
all_fIDR = zeros(1,proteins);
all_ffuzzy = zeros(1,proteins);
all_fresidual = zeros(1,proteins);


figure; hold on
C = get(gca, 'ColorOrder');

for p = 1:proteins
    UniProtID = data{p,1};
    x = str2double(data{p,2});
    all_fIDR(p) = x;
    y = str2double(data{p,3});
    all_ffuzzy(p) = y;
    z = str2double(data{p,4});
    all_fresidual(p) = z;
    psize = str2double(data{p,5});
    nd = str2double(data{p,6});
    xyz = [x,y,1 - z];
    radius = 0.001*psize^(1/3);
    nIFR = nd;
    if nIFR == 0
        rgb = [0,0,0];
    else
        if nIFR > 5
            nIFR = 5;
        end
        rgb = C(nIFR,:);
    end
    plot_pod_sphere(xyz,radius,rgb,UniProtID);
end

axis equal
axis([0,1,0,1,0,1]);
view(35,30);
lighting gouraud 
material shiny

xlabel('f_{IDR}');
ylabel('f_{fuzzy}');
zlabel('1 - f_{residual}');

