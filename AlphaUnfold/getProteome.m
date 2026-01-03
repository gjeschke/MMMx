function UniProtIDs = getProteome(proteomeID,is_reference)
%    UniProtIDs = getProteome(proteomeID,is_reference)
%           Gets UniProt IDs for a proteome
%
% INPUT:
% proteome_id   identifier, such as UP000000625 for the E. coli reference
%               proteome
% is_reference  flag that suppresses query for reviewed UniProt entries and
%               removal of unreviwed entries, defaults to false
%
% OUTPUT:
% UniProtIDs    cell string of UniProt identifiers of the proteins in this
%               proteome

if ~exist('is_reference','var')
    is_reference = false;
end

UniProtIDs = cell(1,100000);
np = 0;

fprintf('Starting download of UniProt IDs...\n');

baseURL = 'https://rest.uniprot.org/uniprotkb/search';
    query = sprintf('proteome:%s', proteomeID);
    fields = 'accession'; % We only want the protein identifiers
    format = 'tsv';
    size = 500; % Maximum allowed per request
    
    % Construct the initial URL
    baseURL = sprintf('%s?query=%s&fields=%s&format=%s&size=%d', ...
                  baseURL, query, fields, format, size);

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
        if startsWith(lines{1}, 'Entry')
            lines = lines(2:end);
        end

        % Add IDs
        if ~isempty(lines)
            fprintf('Got %d IDs\n', numel(lines));
            for i = 1:numel(lines)
                np = np + 1;
                UniProtIDs{np} = lines{i};
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

UniProtIDs = UniProtIDs(1:np);

if ~is_reference
    fprintf(1,'Filtering reviewed entries only\n');
    nr = 0;
    for n = 1:np
        info = getInfoFromUniProt(UniProtIDs{n});
        if contains(info.entryType,'Swiss-Prot')
            nr = nr + 1;
            UniProtIDs{nr} = UniProtIDs{n};
        end
        if mod(n,100) == 0
            fprintf(1,'%i of %i entries processed\n',n,np);
        end
    end
    UniProtIDs = UniProtIDs(1:nr);
else
    nr = np;
end

fprintf(1,'%i proteins in this proteome\n',nr);

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

