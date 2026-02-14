function h = plot_pod_isosurface(data_cube,rgb,n,tag,level)

if ~exist('level','var') || isempty(level)
    level = 0.75;
end
X = linspace(0,1,101);
Y = X;
Z = X;
sdens = sum(sum(sum(data_cube)));
max_population = max(max(max(data_cube)));
isovalue = max_population/100;
for kl = 0.01:0.01:0.99
    mask = (data_cube >= kl*max_population);
    test = sum(sum(sum(mask.*data_cube)));
    if test <= level*sdens
        isovalue = kl*max_population;
        break;
    end
end

h = isosurface(X,Y,Z,data_cube,isovalue);
p = patch(h);
set(p,'FaceColor',rgb);  
set(p,'EdgeColor','none');

p.UserData.tag = tag;
p.UserData.n = n;
p.UserData.selected = false;
p.UserData.rgb = rgb;
p.ButtonDownFcn = @isosurface_clicked;
