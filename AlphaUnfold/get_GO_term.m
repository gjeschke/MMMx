function [GO_term,info] = get_GO_term(go_id)

baseURL = 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/';
url = sprintf('%sGO:%s',baseURL,go_id);

% Call the web API (MATLAB R2016b+)
data = webread(url);

% QuickGO returns a struct; term name is in data.results(1).name
GO_term = data.results(1).name;
info = data.results;