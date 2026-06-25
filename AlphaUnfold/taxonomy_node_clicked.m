function taxonomy_node_clicked(cb,eventdata)
% taxonomy_node_clicked(cb,eventdata)
%
% Function executed when a user clicks on a marker that represents one
% taxonomy node in a taxonomy tree

url = sprintf('https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=%i',cb.UserData.TaxonId);
if isfield(cb.UserData,'proteins')
    proteins = cb.UserData.proteins;
else
    proteins = NaN;
end

switch eventdata.Button
    case 1 % left button
        fprintf(1,'%s (%s) with %i proteins\n',cb.UserData.name,cb.UserData.proteome, proteins);
        title(sprintf('%s with %i proteins, %s = %.3f',cb.UserData.name,proteins,...
        cb.DataTipTemplate.DataTipRows(2).Label,cb.DataTipTemplate.DataTipRows(2).Value));
        fprintf(1,'%s = %.3f\n',cb.DataTipTemplate.DataTipRows(2).Label,cb.DataTipTemplate.DataTipRows(2).Value);
    case 2 % middle button
        web(url, '-browser');
end