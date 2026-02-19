function uniref_id = uniprot2uniref100_search(uniprot_id)
% UNIPROT2UNIREF100_SEARCH  Get UniRef100 cluster ID for any UniProt accession
%
%   uniref_id = uniprot2uniref100_search('P43403')
%
%   Input:
%       uniprot_id : character vector or string (any UniProt accession)
%   Output:
%       uniref_id : character vector with UniRef100 cluster ID, or empty [] if not found
%
%   This function searches the UniRef100 database for clusters containing
%   the given UniProt accession as a member.

    % Clean the input
    uniprot_id = strtrim(char(uniprot_id));
    if isempty(uniprot_id)
        uniref_id = [];
        return;
    end
    
    % Construct the search URL for UniRef100
    % Query: find UniRef100 clusters where uniprot_id is in the members list
    % Format: JSON for easy parsing
    base_url = 'https://rest.uniprot.org/uniref/search';
    query = sprintf('(uniprot_id:%s) AND (identity:1.0)', uniprot_id);
    
    % URL encode the query
    encoded_query = strrep(query, ' ', '%20');
    encoded_query = strrep(encoded_query, ':', '%3A');
    encoded_query = strrep(encoded_query, '(', '%28');
    encoded_query = strrep(encoded_query, ')', '%29');
    
    url = [base_url '?query=' encoded_query '&format=json&size=1'];
    
    try
        % Fetch the search results
        json_str = urlread(url);
        
        % Parse JSON
        results = jsondecode(json_str);
        
        % Check if we got any results
        if isfield(results, 'results') && ~isempty(results.results)
            % Extract the UniRef ID from the first result
            first_result = results.results(1);
            if isfield(first_result, 'id')
                uniref_id = first_result.id;
            else
                uniref_id = [];
            end
        else
            uniref_id = [];
        end
        
    catch ME
        % Check if it's a 404 or other error
        if contains(ME.message, '404') || contains(ME.message, 'Not Found')
            uniref_id = [];
        else
            warning('Error querying %s: %s', uniprot_id, ME.message);
            uniref_id = [];
        end
    end
end