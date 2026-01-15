function protein_clicked(cb,~)
% gobject_clicked(cb,eventdata)
%
% Function executed when a user clicks on a sphere that represents one
% protein in the proteome

url = sprintf('https://www.uniprot.org/uniprotkb/%s/entry',cb.UserData.tag);
web(url, '-browser');
url = sprintf('https://alphafold.ebi.ac.uk/entry/AF-%s-F1',cb.UserData.tag);
web(url, '-browser');
if isfield(cb.UserData,'characteristics') && sum(isnan(cb.UserData.characteristics)) == 0
    fprintf(1,'f(IDR) = %5.3f\n',cb.UserData.characteristics(1));
    fprintf(1,'f(fuzzy) = %5.3f\n',cb.UserData.characteristics(2));
    fprintf(1,'1 - f(residual) = %5.3f\n',cb.UserData.characteristics(3));
end

if isfield(cb.UserData,'GO_ids') && ~isempty(cb.UserData.GO_ids)
    all_ids = split(cb.UserData.GO_ids,';');
    for k = 1:length(all_ids)
        url = sprintf('http://purl.obolibrary.org/obo/GO_%s',all_ids{k});
        web(url, '-browser');
    end
end
