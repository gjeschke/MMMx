function np = getUniProtIDs(taxonomy_id,filename,subclass)
    % Gets UniProt IDs for Riboviria using urlread (simpler, returns text)
    % INPUT: filename (string) - e.g., 'riboviria_ids.txt'

    if ~exist('subclass','var')
        subclass = '';
    end

    if ~exist('filename','var') || isempty(filename)
        filename = 'UniProt_ids.dat';
    end

    if isempty(subclass)
        fid = fopen(filename,'w');
    else
        fid = fopen(filename,'a');
    end

    if fid == -1
        fprintf('Error opening file.\n');
        return;
    end
    
    np = 0;

    fprintf('Starting download of UniProt IDs...\n');
    
    % Base URL
    baseURL = sprintf('https://rest.uniprot.org/uniprotkb/search?format=list&query=(taxonomy_id:%i)&size=500',taxonomy_id);
    
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
                for i = 1:numel(lines)
                    np = np + 1;
                    if isempty(subclass)
                        fprintf(fid, '%s\n', lines{i});
                    else
                        fprintf(fid, '%s,%s\n', lines{i},subclass);
                    end
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
    
    % Save results
    fclose(fid);
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

