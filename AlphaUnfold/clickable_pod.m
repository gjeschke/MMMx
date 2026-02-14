function my_light = clickable_pod(options)

if ~exist('options','var') || ~isfield(options,'opaqueness')
    options.opaqueness = [];
end

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
all_psize = zeros(1,proteins);


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
    all_psize(p) = psize; 
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
    if isempty(options.opaqueness)
        opaqueness = 1;
    else
        opaqueness = str2double(data{p,options.opaqueness});
    end
    plot_pod_sphere(xyz,radius,rgb,UniProtID,opaqueness,[all_fIDR(p),all_ffuzzy(p),1 - all_fresidual(p)],data{p,7},psize);
end

axis equal
grid on;
axis([0,1,0,1,0,1]);
view(35,30);
lighting gouraud 
material shiny
my_light = camlight;

xlabel('f_{IDR}');
ylabel('f_{fuzzy}');
zlabel('1 - f_{conditional}');

fprintf(1,'Mean f(IDR): %6.3f\n',sum(all_psize.*all_fIDR)/sum(all_psize));
fprintf(1,'Mean f(fuzzy): %6.3f\n',sum(all_psize.*all_ffuzzy)/sum(all_psize));
fprintf(1,'Mean f(residual): %6.3f\n',sum(all_psize.*all_fresidual)/sum(all_psize));

