function isosurface_clicked(cb,~)
% GO_clicked(cb,eventdata)
%
% Function executed when a user clicks on a sphere that represents one
% gene ontology


if ~cb.UserData.selected
    fprintf(1,'Selected: %s with %i proteins\n',cb.UserData.tag,cb.n);

    rgb = [0.9290    0.6940    0.1250];
    set(cb,'FaceColor',rgb);  
else

    set(cb,'FaceColor',cb.UserData.rgb);
end

