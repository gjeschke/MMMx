function upids = get_proteome_by_taxonomy_id(taxonomy_id)

upids = {};

url = sprintf('https://rest.uniprot.org/proteomes/search?query=%i',taxonomy_id);
[response, status] = urlread(url); %#ok<URLRD> 
if status ~= 1 || isempty(response)
    return
end
response = jsondecode(response);
proteome_id = '';
proteins = 0;
is_reference = false;
if length(response.results) == 1
    if iscell(response.results.taxonLineage)
        this_id = response.results.taxonLineage{end}.taxonId;
    else
        this_id = response.results.taxonLineage(end).taxonId;
    end
    if this_id == taxonomy_id % because this server also serves unrelated stuff
        proteome_id = response.results.id;
    else
        return
    end
else
    for k = 1:length(response.results)
        if iscell(response.results)
            this_response = response.results{k};
        else
            this_response = response.results(k);
        end
        if iscell(this_response.taxonLineage)
            this_id = this_response.taxonLineage{end}.taxonId;
        else
            this_id = this_response.taxonLineage(end).taxonId;
        end
        if this_id == taxonomy_id % because this server also serves unrelated stuff
            if ~is_reference % if we do not yet have a reference proteome, use this one
                if this_response.proteinCount > proteins % but only if it has more proteins than an already stored proteome
                    proteins = this_response.proteinCount;
                    proteome_id = this_response.id;
                end
            end
            ptype = this_response.proteomeType;
            if contains(lower(ptype),'reference') || contains(lower(ptype),'representative')
                if this_response.proteinCount > proteins
                    proteins = this_response.proteinCount;
                    proteome_id = this_response.id;
                    is_reference = true;
                end
            end
        end
    end
end
if isempty(proteome_id)
    fprintf(2,'No proteome found for taxonomy ID %i\n',taxonomy_id);
else
    upids = getProteome(proteome_id,true);
end


