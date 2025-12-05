function protein_clicked(cb,~)
% gobject_clicked(cb,eventdata)
%
% Function executed when a user clicks on a sphere that represents one
% protein in the proteome

url = sprintf('https://www.uniprot.org/uniprotkb/%s/entry',cb.UserData.tag);
web(url, '-browser');
url = sprintf('https://alphafold.ebi.ac.uk/entry/AF-%s-F1',cb.UserData.tag);
web(url, '-browser');
