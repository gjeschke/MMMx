function plot_pod_octahedron(xyz,radius,rgb,tag,alpha,characteristics,GO_ids,residues)

if ~exist('tag','var')
    tag = 'tag';
end

if ~exist('characteristics','var')
    characteristics = nan(1,3);
end

if ~exist('alpha','var')
    alpha = 1;
end

if ~exist('GO_ids','var')
    GO_ids = '';
end


obj = plot_octahedron(xyz, 2.5*radius, 2.5*radius, 2.5*radius, rgb);

set(obj, 'FaceColor', rgb, 'EdgeColor', 'none', 'FaceAlpha',alpha,'FaceLighting','gouraud','Clipping','off');
set(obj, 'CDataMapping','direct','AlphaDataMapping','none');
obj.UserData.tag = tag;
obj.UserData.characteristics = characteristics;
obj.UserData.GO_ids = GO_ids;
obj.UserData.residues = residues;
obj.ButtonDownFcn = @protein_clicked;

end
            