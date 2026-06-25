function GO_node_clicked(cb,eventdata)
% GO_node_clicked(cb,eventdata)
%
% Function executed when a user clicks on a marker that represents one
% GO node in a gene ontology tree

url = sprintf('http://purl.obolibrary.org/obo/GO_%s',cb.UserData.GO_id);
if isfield(cb.UserData,'proteins')
    proteins = cb.UserData.proteins;
else
    proteins = NaN;
end

switch eventdata.Button
    case 1 % left button
        web(url, '-browser');
    case 2 % middle button
        fprintf(1,'%s with %i proteins\n',cb.UserData.name,proteins);
        fprintf(1,'%s = %.3f\n',cb.DataTipTemplate.DataTipRows(2).Label,cb.DataTipTemplate.DataTipRows(2).Value);
end