function characteristics = pod_parameters(UPID_file,options)
% POD_PARAMETER(UPID_file,options)
% Determines some disorder characteristics for a list of UniProt
% identifiers and writes a file with the characteristics for individual
% proteins as well as files for all gene ontology identifiers
%
% Input:
% UPID_file     .csv file with UniProt identifiers for the whole proteome,
%               UniProt ID must be in the first column
% options       run options, struct with fields
%               .threads    number of parallel threads, defaults to 50
%               .input      input description for avoiding recomputation or
%                           repeated database query, 0 means that a value
%                           needs to be recomputed, defaults to all 0, the
%                           numbers denote columns in the .csv input file,
%                           -1 means that the item is not present in input,
%                           but should not be computed, but left empty in
%                           output
%                           .fIDR   fraction of residues in IDRs
%                           .ffuzzy fraction of moderately certain residues
%                                   in IFR
%                           .fcondi fraction of residues in IDRs that
%                                   probably fold conditionally
%                           .size   length of the sequence
%                           .go_ids gene ontology identifiers
%                           .taxon  taxonomy (lineage) of organism
%                           .dynam  dynamics predictor derived from PAE
%               .go_distributed     distribute proteins into files for
%                                   individual gene ontology identifiers,
%                                   flag, defaults to true, prolongs
%                                   computation time significantly for
%                                   large datasets
%
% Output:
% characteristics   disorder characteristics, struct with fields
%                   .fIDR                   fraction of residues in IDRs
%                   .fuzziness              fuzziness of IFRs
%                   .conditional_structure  conditional structure of IDRs
%                   .sequence_lengths       histogram of sequence lengths
%                   .proteins               number of proteins
%                   .oversized              number of proteins with more
%                                           than 5000 residues, not 
%                                           contained in .sequence_lengths
%
% G. Jeschke, 2026


ksi = 0.75; % steepnees parameter for sigmoid function for dynamics prediction 
paeh = 30; % half saturation for sigmoid function for dynamics prediction

if ~exist('options','var') || ~isfield(options,'threads')
    options.threads = 55;
end

if ~isfield(options,'input') || isempty(options.input)
    options.input.fIDR = 0;
    options.input.fuzzy = 0;
    options.input.fcondi = 0;
    options.input.size = 0;
    options.input.go_ids = 0;
    options.input.taxon = 0;
    options.input.dynam = 0;
end

if ~isfield(options,'go_distributed') || isempty(options.go_distributed)
    options.go_distributed = true;
end

[path,fname,ext] = fileparts(UPID_file);
if isempty(ext)
    ext = '.dat';
end
iname = fullfile(path,strcat(fname,ext));
oname = fullfile(path,sprintf('%s_pod_characteristics.csv',fname));

infopoint = 50000;
uniprot_ids = cell(options.threads,1);

all = zeros(1,101);
fuzziness = zeros(1,101);
conditional_structure = zeros(1,101);
dynamics = zeros(1,101);
all_sizes = zeros(1,5000);
oversized = zeros(1,50000);
n_oversized = 0;
fid = fopen(iname);
ofid = fopen(oname,'w');

proteins = 0;
AF_options.structure = false;
AF_options.pLDDT = false;
AF_options.pae = true;

cap = 31;

domain_options.threshold = 0.65*cap;
domain_options.unify = 15;
domain_options.minsize = 15;
domain_options.minlink = 3;
domain_options.local = 15;
domain_options.interact = 1;


orphans = fopen('orphans.csv','w');

k = 0;
AF_failed = 0;
tic,
tline0 = '';
while 1
    for t = 1:options.threads
        tline = fgetl(fid);
        if strcmp(tline,tline0)
            % uniprot_ids{t} = 'xxx';
            % continue;
        else
            tline0 = tline;
        end
        if ~ischar(tline)
            options.threads = t-1;
            break
        end
        args = split(tline,',');
        uniprot_ids{t} = args{1}; 
        if options.input.fIDR > 0
            fIDR = str2double(args{options.input.fIDR});
        elseif options.input.fIDR < 0
            fIDR = [];
        end
        if options.input.fuzzy > 0
            ffuzzy = str2double(args{options.input.fuzzy});
        elseif options.input.fuzzy < 0
            ffuzzy = [];
        end
        if options.input.fcondi > 0
            fcondi = str2double(args{options.input.fcondi});
        elseif options.input.fcondi < 0
            fcondi = [];
        end
        if options.input.size > 0
            psize = str2double(args{options.input.psize});
        elseif options.input.size < 0
            psize = [];
        end
        if options.input.go_ids > 0
            go_ids = args{options.input.go_ids};
        elseif options.input.go_ids < 0
            go_ids = '';
        end
        if options.input.taxon > 0
            taxonomy = args{options.input.taxon};
        elseif options.input.taxon < 0
            taxonomy = '';
        end
        if options.input.dynam > 0
            fdynam = str2double(args{options.input.fdynam});
            AF_options.pae = false;
        elseif options.input.dynam < 0
            fdynam = [];
            AF_options.pae = false;
        end
    end
    if isempty(options.threads) 
        break
    end
    datasets = cell(options.threads,1);
    query_part = strjoin(strcat('accession:', uniprot_ids), ' OR ');
    encoded_query = urlencode(query_part);
    url = sprintf('https://rest.uniprot.org/uniprotkb/search?query=%s&fields=organism_name,go_id&format=json&size=500', encoded_query);
    results = urlread(url); %#ok<URLRD>
    info = jsondecode(results);
    % info = request(url);
    if isempty(info)
        continue
    end
    uniprot_info = info.results;
    parfor t = 1:options.threads %#ok<PFUIXW> 
        try
            entity = get_AF(uniprot_ids{t},AF_options);
            datasets{t} = entity;
        catch
            datasets{t} = [];
        end
    end
    for u = 1:length(uniprot_info)
        if iscell(uniprot_info)
            info = uniprot_info{u};
        else
            info = uniprot_info(u);
        end
        upid = info.primaryAccession;
        t = 0;
        for tu = 1:options.threads
            if strcmpi(upid,uniprot_ids{tu})
                t = tu;
                break
            end
        end
        if t == 0
            break
        end
        entity = datasets{t};
        if isempty(entity)
           AF_failed = AF_failed + 1;
           continue
        end
        AF_info = entity.AF_info;
        if isempty(AF_info)
           AF_failed = AF_failed + 1; 
        end
        if length(AF_info) > 1
            AF_info = AF_info(1);
        end
        if isempty(info) || isempty(AF_info)
            continue
        end
        if iscell(AF_info)
            AF_info = AF_info{1};
        end
        if ~options.input.dynam
            pae = entity.pae;
            domains = get_domains(pae,domain_options);
            [nd,~] = size(domains);
            pae = atan(ksi*(pae-paeh))/pi + 1/2;
            offset = atan(-ksi*paeh)/pi + 1/2;
            pae = pae - offset;
            sc = atan(ksi*(2*cap/3-paeh))/pi + 1/2 - offset;
            pae = pae/sc;
            mean_PAE = 0;
            nsum = 0;
            for d = 1:nd
                mean_PAE = mean_PAE + (1 + domains(d,2) - domains(d,1))*mean(pae(domains(d,1):domains(d,2),domains(d,1):domains(d,2)),'all');
                nsum = nsum + 1 + domains(d,2) - domains(d,1);
            end
            fdynam = mean_PAE/nsum;
        end
        if ~isfield(info,'organism')
            continue
        end
        if length(info.organism.lineage) < 5
            fprintf(orphans,'%s,%s\n',uniprot_ids{t},strjoin(info.organism.lineage, ','));
            continue
        end
        proteins = proteins + 1;
        if ~options.input.size
            psize = length(AF_info.sequence);
        end
        if psize > length(all_sizes)
            n_oversized = n_oversized + 1;
            oversized(n_oversized) = psize;
        else 
            all_sizes(psize) = all_sizes(psize) + 1;
        end
        if options.input.fIDR == 0 || options.input.ffuzzy == 0 || options.input.fcondi == 0
            fIDR = AF_info.fractionPlddtVeryLow + AF_info.fractionPlddtLow;
            index = 1 + round(100*fIDR);
            all(index) = all(index) + 1;
            if floor(fIDR*psize) > 15 % at least 15 disordered residues
                fcondi = AF_info.fractionPlddtLow/(AF_info.fractionPlddtVeryLow + AF_info.fractionPlddtLow);
                index = 1 + round(100*fcondi);
                conditional_structure(index) = conditional_structure(index) + 1;
            else
                fcondi = 1;
            end
            if floor((1-fIDR)*psize) > 15 % at least 15 residues in IFRs
                ffuzzy = AF_info.fractionPlddtConfident/(AF_info.fractionPlddtVeryHigh + AF_info.fractionPlddtConfident);
                index = 1 + round(100*ffuzzy);
                fuzziness(index) = fuzziness(index) + 1;
            else
                ffuzzy = 0;
            end
        end
        if ~options.input.taxon
            taxonomy = strjoin(info.organism.lineage, ';');
        end
        if ~options.input.go_ids
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
                        if isempty(go_ids)
                            go_ids = go_info{2};
                        else
                            go_ids = sprintf('%s;%s',go_ids,go_info{2});
                        end
                        if options.go_distributed
                            goname = sprintf('GO_%s.csv',go_info{2});
                            gfid = fopen(goname,'a');
                            fprintf(gfid,'%s,%5.3f,%6.3f,%6.3f,%i,%s,%s,%6.3f\n',...
                                uniprot_ids{t},fIDR,ffuzzy,fcondi,psize,go_info{2},taxonomy,fdynam);
                            fclose(gfid);
                        end
                    end
                end
            end
        end
        fprintf(ofid,'%s,%5.3f,%6.3f,%6.3f,%i,%s,%s,%6.3f\n',...
            uniprot_ids{t},fIDR,ffuzzy,fcondi,psize,go_ids,taxonomy,fdynam);
    end
    if t == 0
        break
    end
    k = k + options.threads;
    if mod(k,infopoint) == 0
        fprintf(1,'%6.3f%% of proteins read in %i attempts\n',100*proteins/k,k);
    end
    if ~ischar(tline)
        break
    end
end
toc,
oversized = oversized(1:n_oversized);

fclose(fid);
fclose(orphans);
fclose(ofid);

characteristics.fIDR = all;
characteristics.fuzziness = fuzziness;
characteristics.conditional_structure = conditional_structure;
characteristics.sequence_lengths = all_sizes;
characteristics.proteins = proteins;
characteristics.oversized = oversized;
characteristics.dynamics = dynamics;

fprintf(1,'%i proteins were read\n',proteins);
fprintf(1,'For %i proteins, AlphaFold information could not be read\n',AF_failed);

% function data = request(url)
% 
% import matlab.net.*
% import matlab.net.http.*
% 
% persistent options 
% 
% data = [];
% 
% % import matlab.net.*
% % import matlab.net.http.*
% 
% % 1. Create the Request Message
% request = RequestMessage();
% 
% % 2. Define HTTP Options (Crucial for performance)
% if isempty(options)
%     options = HTTPOptions(...
%         'ConnectTimeout', 10, ...
%         'KeepAliveTimeout', Inf, ... % Reuse connection
%         'ConvertResponse', false ... % Skip auto-conversion
%         );
% end
% % 3. Create the HTTP Client correctly
% % client = matlab.net.http.HttpClient();
% 
% % 4. Send the request
% % response = client.send(request, url, options);
% response = request.send(url, options);
% 
% % 5. Manually handle the response (e.g., JSON)
% if ~isempty(response.Body.Data)
%     rawData = native2unicode(response.Body.Data');
%     try
%         data = jsondecode(rawData);
%     catch
%         data = [];
%     end
% end