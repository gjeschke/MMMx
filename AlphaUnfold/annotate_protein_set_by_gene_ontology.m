function annotation = annotate_protein_set_by_gene_ontology(UPID_file)
% Distributes UniProt identifiers into separate files named by the
% corresponding gene ontology identifiers, one protein can have several GO
% identifiers
%
% Input:
% UPID_file     file with UniProt identifiers for the whole proteome, must
%               be a .csv file and the UniProt identifier must be in the
%               first column
%
% G. Jeschke, 2026

annotation(250000).upid = '';
annotation(250000).goid = 0;
annotation(250000).taxonid = 0;
annotation(250000).hostid = 0;

annotated = 0;

batch_size = 50;

[path,fname,ext] = fileparts(UPID_file);
if isempty(ext)
    ext = '.csv';
end
iname = fullfile(path,strcat(fname,ext));
oname = fullfile(path,sprintf('%s_with_GO.csv',fname));
ofid = fopen(oname,'w');

infopoint = 1000;
uniprot_ids = cell(batch_size,1);
taxonomy_ids = zeros(batch_size,1);
host_ids = zeros(batch_size,1);

fid = fopen(iname);

proteins = 0;

batch = 0;

while 1
    for t = 1:batch_size
        tline = fgetl(fid);
        if ~ischar(tline)
            batch_size = t - 1;
            break
        end
        args = split(tline,',');
        uniprot_ids{t} = args{1}; 
        taxonomy_ids(t) = str2double(args{2}); 
        host_ids(t) = str2double(args{3}); 
    end
    if batch_size > 0
        batch = batch + 1;
    else
        break
    end
    uniprot_ids = uniprot_ids(1:batch_size);
    query_part = strjoin(strcat('accession:', uniprot_ids), ' OR ');
    encoded_query = urlencode(query_part);
    url = sprintf('https://rest.uniprot.org/uniprotkb/search?query=%s&fields=organism_name,go_id&format=json&size=500', encoded_query);
    data = urlread(url); %#ok<URLRD> 
    info = jsondecode(data);
    %info = request(url);
    if isempty(info)
        fprintf(2,'Warning no data in batch %i (after %i proteins)\n',batch,proteins);
        continue
    end
    uniprot_info = info.results;
    for u = 1:length(uniprot_info)
        info = uniprot_info(u);
        upid = info.primaryAccession;
        t = 0;
        for tu = 1:batch_size
            if strcmpi(upid,uniprot_ids{tu})
                t = tu;
                break
            end
        end
        if t == 0
            fprintf(2,'UniProt ID %s is surplus information\n',upid);
            break
        end
        if ~isfield(info,'organism')
            fprintf(2,'Organism information missing for UniProt ID %s\n',upid)
            continue
        end
        proteins = proteins + 1;
        taxonomy = strjoin(info.organism.lineage, ';');
        go_ids = '';
        if isfield(info,'uniProtKBCrossReferences')
            for kg = 1:length(info.uniProtKBCrossReferences)
                clear('go_info');
                if iscell(info.uniProtKBCrossReferences)
                    if strcmpi(info.uniProtKBCrossReferences{kg}.database,'GO')
                        go_info = split(info.uniProtKBCrossReferences{kg}.id,':');   
                    end
                else
                    if strcmpi(info.uniProtKBCrossReferences(kg).database,'GO')
                        go_info = split(info.uniProtKBCrossReferences(kg).id,':');   
                    end
                end
                if exist('go_info','var')
                    goname = sprintf('GO_%s.csv',go_info{2});
                    go_ids = sprintf('%s;%s',go_ids,go_info{2});
                    gfid = fopen(goname,'a');
                    fprintf(gfid,'%s,%s,%s\n',...
                        uniprot_ids{t},go_info{2},taxonomy);
                    fclose(gfid);
                    annotated = annotated + 1;
                    annotation(annotated).upid = uniprot_ids{t};
                    annotation(annotated).goid = str2double(go_info{2});
                    annotation(annotated).taxonid = taxonomy_ids(t);
                    annotation(annotated).hostid = host_ids(t);                
                end
            end
        end
        if ~isempty(go_ids) % remove leading semicolon
            go_ids = go_ids(2:end);
        end
        fprintf(ofid,'%s,%s,%s\n',...
            uniprot_ids{t},go_ids,taxonomy);
    end
    if t == 0
        break
    end
    if mod(proteins,infopoint) == 0
        fprintf(1,'%i proteins were read\n',proteins);
    end
    if ~ischar(tline)
        break
    end
end

fclose(fid);
fclose(ofid);

fprintf(1,'Finished. %i proteins were read\n',proteins);
fprintf(1,'%i annotations were made\n',annotated);

annotation = annotation(1:annotated);

