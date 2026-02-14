function [UniProtIDs,organism] = getUniProtIDs_by_TaxID(taxID)
    % Gets UniProt IDs for Riboviria using urlread (simpler, returns text)
    % INPUT: filename (string) - e.g., 'riboviria_ids.txt'

    baseURL = sprintf('https://rest.uniprot.org/taxonomy/search?query=%i&format=json',taxID);
    response = urlread(baseURL); %#ok<URLRD> 
    response = jsondecode(response);
    organism = response.results.scientificName;

    UniProtIDs = cell(100000,1);
    proteins = 0;

    fprintf('Starting download of UniProt IDs...\n');
    
    % Base URL
    baseURL = sprintf('https://rest.uniprot.org/uniprotkb/search?query=organism_id:%i&format=tsv&fields=accession&size=500',taxID);
    
    nextURL = baseURL;  % Start with base URL
    batchCount = 0;
    
    while ~isempty(nextURL)
        batchCount = batchCount + 1;
        fprintf('Batch %d... ', batchCount);
        
        try
            % 1. Get data using urlread (returns text)
             [response, status] = urlread(nextURL); %#ok<URLRD>
            
            if status ~= 1 || isempty(response)
                fprintf('Failed to get data\n');
                break;
            end
            
            % 2. Process the text response
            lines = strsplit(response, {'\r\n', '\n', '\r'});
            lines = lines(~cellfun('isempty', lines));
            
            % Remove header line if present
            if batchCount == 1 && startsWith(lines{1}, 'accession')
                lines = lines(2:end);
            end
            
            % Add IDs
            if ~isempty(lines)
                fprintf('Got %d IDs\n', numel(lines));
                for i = 2:numel(lines)
                    proteins = proteins + 1;
                    UniProtIDs{proteins} = lines{i};
                end
            else
                fprintf('No IDs in batch\n');
            end
            
            % 3. Check for more data using Java
            nextURL = checkForMoreJava(nextURL);
            
            % Small pause
            pause(0.3);
            
        catch ME
            fprintf('Error: %s\n', ME.message);
            fprintf('Saving what we have...\n');
            break;
        end
    end

    UniProtIDs = UniProtIDs(1:proteins);
    fprintf(1,'Retrieved %i UniProt identifiers for taxonomy ID %i\n',proteins,taxID);
    
end

function nextURL = checkForMoreJava(currentURL)
    % Check if there are more results using Java
    import java.net.* java.io.*
    
    nextURL = '';
    
    try
        
        % Make a HEAD request
        urlObj = URL(currentURL);
        conn = urlObj.openConnection();
        conn.setRequestMethod('HEAD');
        conn.setConnectTimeout(3000);
        conn.setReadTimeout(3000);
        
        % Get Link header
        linkHeader = conn.getHeaderField('Link');
        
        if ~isempty(linkHeader)
            % Look for next page
            pattern = '<([^>]+)>;\s*rel="next"';
            tokens = regexp(char(linkHeader), pattern, 'tokens');
            
            if ~isempty(tokens)
                nextURL = tokens{1}{1};
            end
        end
    catch
        % No more pages
    end
end

