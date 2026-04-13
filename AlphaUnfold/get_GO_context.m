function [GO_term,info] = get_GO_context(go_id)

context = '/ancestors?relations=is_a%2Cpart_of%2Coccurs_in%2Cregulates';

baseURL = 'https://www.ebi.ac.uk/QuickGO/services/ontology/go/terms/';
url = sprintf('%sGO:%s%s',baseURL,go_id,context);

% Call the web API (MATLAB R2016b+)
data = webread(url);

% QuickGO returns a struct; term name is in data.results(1).name
GO_term = data.results(1).name;
info = data.results;