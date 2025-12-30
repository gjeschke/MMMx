function info = getInfoFromUniProt(uniprot_id)
% GETINFOFROMUNIPROT Fetches organism name and lineage for a single UniProt ID.
%   INFO = GETINFOFROMUNIPROT(UNIPROT_ID) returns the
%   full entry information, or an empty string on error.
%

    % Initialize the output
    info = '';
    
    % Construct the URL (as in your working code)
    url = sprintf('https://rest.uniprot.org/uniprotkb/search?query=accession:%s&fields=organism_name&format=json', uniprot_id);
    
    try
        % 1. Fetch the raw JSON text
        rawText = urlread(url); %#ok<URLRD>
        
        % 2. Parse the JSON
        jsonData = jsondecode(rawText);
        
        % 3. Extract the name, checking the structure exists
        if isfield(jsonData, 'results') && ~isempty(jsonData.results)
            % The field is nested: results.organism.scientificName
            info = jsonData.results;
        else
            warning('UniProt ID "%s" not found or has no entry data.', uniprot_id);
        end
        
    catch ME
        % Handle specific common errors
        switch true
            case contains(ME.message, '404')
                warning('UniProt ID "%s" does not exist (404 Not Found).', uniprot_id);
            case contains(ME.message, 'timed out') || contains(ME.message, 'Timeout')
                warning('Request for "%s" timed out. Server may be busy.', uniprot_id);
            otherwise
                warning('Failed to process "%s". Error: %s', uniprot_id, ME.message);
        end
    end
end