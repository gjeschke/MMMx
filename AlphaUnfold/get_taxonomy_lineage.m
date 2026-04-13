function [lineage,organism] = get_taxonomy_lineage(taxon_id)
% lineage = get_taxonomy_lineage(taxon_id)
%
% Input:
% taxon_id      NCBI taxonomy identifier
%
% Output:
% lineage   taxonomy lineage in terms of identifiers starting from root (1)
%           ends at parent, empty if the query to NCBI fails
%
% G. Jeschke, 4.4.2026

options = weboptions("ContentType","json","Timeout",30);
url = sprintf('https://api.ncbi.nlm.nih.gov/datasets/v2/taxonomy/taxon/%i',taxon_id);
trials = 10;
result = '';
while trials > 0
    try
        result = webread(url,options);
    catch
        trials = trials - 1;
    end
    if exist('result','var') && ~isempty(result)
        break
    end
end
try
    lineage = result.taxonomy_nodes.taxonomy.lineage;
catch
    lineage = [];
end

try
    organism = result.taxonomy_nodes.taxonomy.organism_name;
catch
    organism = '';
end

