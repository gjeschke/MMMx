function plot_proteomes(proteome_list)

data = get_csv(proteome_list,'data.cell');
[proteins,~] = size(data);
all_fIDR = zeros(1,proteins);
all_ffuzzy = zeros(1,proteins);
all_fresidual = zeros(1,proteins);


figure; hold on
C = get(gca, 'ColorOrder');

for p = 1:proteins
    Organism = data{p,1};
    x = str2double(data{p,3});
    all_fIDR(p) = x;
    y = str2double(data{p,4});
    all_ffuzzy(p) = y;
    z = str2double(data{p,5});
    all_fresidual(p) = z;
    psize = str2double(data{p,2});
    ptype = str2double(data{p,6});
    xyz = [x,y,1 - z];
    radius = 0.02; % 0.001*psize^(1/3);
    if ptype > 5
        rgb = C(ptype-5,:);
    else
        rgb = C(ptype,:);
    end
    % if ptype == 0
    %     rgb = [0,0,0];
    % else
    %     if ptype > 5
    %         ptype = 5;
    %     end
    %     rgb = C(ptype,:);
    % end
    if ptype <= 5
        plot_pod_sphere(xyz,radius,rgb,Organism);
    else
        plot_octahedron(xyz,2.5*radius,2.5*radius,2.5*radius,rgb);
    end
end

axis equal
grid on;
axis([0,1,0,1,0,1]);
view(45,20);
lighting gouraud 
material shiny

xlabel('\langle f_{IDR} \rangle');
ylabel('\langle f_{fuzzy} \rangle');
zlabel('\langle 1 - f_{residual} \rangle');
