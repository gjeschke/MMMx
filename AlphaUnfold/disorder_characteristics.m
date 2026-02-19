function characteristics = disorder_characteristics(UPID_file,options)
% DISORDER_CHARACTERISTICS(UPID_file)
% Determines some disorder characteristics for the proteome of an organism
% and writes a file with the characteristics for individual proteins
%
% Input:
% UPID_file     file with UniProt identifiers for the whole proteome
% options       run options, struct with fields
%               .threads    number of parallel threads, defaults to 50
%               .path       path to AF3 predictions, defaults to current
%                           directory
%               .jpath      path to viro3D json scores files, defaults to
%                           current directory
%               .lineage    if true, the entries are distributed into a
%                           directory tree according to lineage, with the
%                           file names derived from species name, defaults
%                           to true
%               .reviewed   if true, only reviewed entries are considered, 
%                           defaults to false 
%
% Output:
% characteristics   disorder characteristics, struct with fields
%                   .fIDR                   fraction of residues in IDRs
%                   .fuzziness              fuzziness of IFRs
%                   .residual_structure     residual structure of IDRs
%                   .sequence_lengths       histogram of sequence lengths
%                   .nIFR                   histogram of the number of IFRs
%                   .proteins               number of proteins
%                   .oversized              number of proteins with more
%                                           than 5000 residues, not 
%                                           contained in .sequence_lengths
%
% G. Jeschke, 2025


if ~exist('options','var') || ~isfield(options,'threads')
    options.threads = 50;
end

if ~isfield(options,'path')
    options.path = '';
end

if ~isfield(options,'jpath')
    options.jpath = '';
end

if ~isfield(options,'lineage') || isempty(options.lineage)
    options.lineage = false;
end

if ~isfield(options,'reviewed') || isempty(options.reviewed)
    options.reviewed = false;
end

[path,fname,ext] = fileparts(UPID_file);
if isempty(ext)
    ext = '.dat';
end
iname = fullfile(path,strcat(fname,ext));
oname = fullfile(path,sprintf('%s_disorder_characteristics.csv',fname));

infopoint = 1000;
uniprot_ids = cell(options.threads,1);
uniprot_info = cell(options.threads,1);
viro3d_ids = cell(options.threads,1);

all = zeros(1,101);
fuzziness = zeros(1,101);
residual_structure = zeros(1,101);
all_sizes = zeros(1,5000);
nIFR = zeros(1,100);
oversized = zeros(1,50000);
n_oversized = 0;
fid = fopen(iname);
if ~options.lineage
    ofid = fopen(oname,'w');
end

proteins = 0;
options.structure = false;
options.pLDDT = true;

orphans = fopen('orphans.csv','w');

k = 0;
tic,
tline0 = '';
while 1
    for t = 1:options.threads
        tline = fgetl(fid);
        if strcmp(tline,tline0)
            uniprot_ids{t} = 'xxx';
            continue;
        else
            tline0 = tline;
        end
        if ~ischar(tline)
            options.threads = t-1;
            break
        end
        args = split(tline,',');
        uniprot_ids{t} = args{1}; 
        if length(args) >= 2
            viro3d_ids{t} = args{2};
        else
            viro3d_ids{t} = 'xxx';
        end
    end
    for t = 1:options.threads % ### parfor
        if ~strcmpi(uniprot_ids{t},'xxx')
            uniprot_info{t} = getInfoFromUniProt(uniprot_ids{t});
        else
            uniprot_info{t} = '';
        end
        ename = sprintf('%s.mat',uniprot_ids{t});
        if ~isempty(options.path) %#ok<PFBNS> 
            ename = fullfile(options.path,ename);
        end
        jname = sprintf('%s*.json',viro3d_ids{t});
        jname = fullfile(options.jpath,jname);
        my_list = dir(jname);
        if exist(ename,'file')
            data = load(ename);
            dataset = data.entity;
            dataset.pae = data.pae;
            dataset.pLDDT = data.pLDDT;
            dataset.AF_info.fractionPlddtVeryLow = sum(data.pLDDT < 50)/length(data.pLDDT);
            dataset.AF_info.fractionPlddtLow = sum((data.pLDDT >= 50) & (data.pLDDT < 70))/length(data.pLDDT);
            dataset.AF_info.fractionPlddtConfident = sum((data.pLDDT >= 70) & (data.pLDDT < 90))/length(data.pLDDT);
            dataset.AF_info.fractionPlddtVeryHigh = sum(data.pLDDT >= 90)/length(data.pLDDT);
            datasets{t} = dataset;
        elseif length(my_list) == 1
            data = jsondecode(fileread(fullfile(my_list.folder,my_list.name)));
            dataset.pae = data.pae;
            dataset.pLDDT = data.plddt;
            dataset.AF_info.fractionPlddtVeryLow = sum(data.plddt < 50)/length(data.plddt);
            dataset.AF_info.fractionPlddtLow = sum((data.plddt >= 50) & (data.plddt < 70))/length(data.plddt);
            dataset.AF_info.fractionPlddtConfident = sum((data.plddt >= 70) & (data.plddt < 90))/length(data.plddt);
            dataset.AF_info.fractionPlddtVeryHigh = sum(data.plddt >= 90)/length(data.plddt);
            datasets{t} = dataset;
        else
            datasets{t} = get_AF(uniprot_ids{t},options);
        end
    end
    for t = 1:options.threads
        info = uniprot_info{t};
        entity = datasets{t};
        if isempty(entity) || isempty(entity.AF_info) || isempty(entity.pae)
            continue
        end
        if ~isfield(info,'organism')
            continue
        end
        if length(info.organism.lineage) < 5
            fprintf(orphans,'%s,%s\n',uniprot_ids{t},strjoin(info.organism.lineage, ','));
            continue
        end
        if options.reviewed && ~contains(info.entryType,'Swiss-Prot')
            continue
        end
        proteins = proteins + 1;
        psize = length(entity.pLDDT);
        if psize > length(all_sizes)
            n_oversized = n_oversized + 1;
            oversized(n_oversized) = psize;
        else 
            all_sizes(psize) = all_sizes(psize) + 1;
        end
        fraction = entity.AF_info.fractionPlddtVeryLow + entity.AF_info.fractionPlddtLow;
        index = 1 + round(100*fraction);
        all(index) = all(index) + 1;
        if floor(fraction*psize) > 15 % at least 15 disordered residues
            fresidual = entity.AF_info.fractionPlddtLow/(entity.AF_info.fractionPlddtVeryLow + entity.AF_info.fractionPlddtLow);
            index = 1 + round(100*fresidual);
            residual_structure(index) = residual_structure(index) + 1;
        else
            fresidual = 1;
        end
        if floor((1-fraction)*psize) > 15 % at least 15 residues in IFRs
            ffuzzy = entity.AF_info.fractionPlddtConfident/(entity.AF_info.fractionPlddtVeryHigh + entity.AF_info.fractionPlddtConfident);
            index = 1 + round(100*ffuzzy);
            fuzziness(index) = fuzziness(index) + 1;
        else
            ffuzzy = 0;
        end
        domains = get_domains(entity.pae);
        [nd,~] = size(domains);
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
                    fprintf(gfid,'%s,%5.3f,%6.3f,%6.3f,%i,%i,%s,%s\n',...
            uniprot_ids{t},fraction,ffuzzy,fresidual,psize,nd,go_info{2},taxonomy);
                    fclose(gfid);
                end
            end
        end
        if ~isempty(go_ids) % remove leading semicolon
            go_ids = go_ids(2:end);
        end
        if options.lineage
            go_home = pwd;
            for tax = 1:4
                tax_term = unblank(info.organism.lineage{tax});
                if ~exist(tax_term,'dir')
                    mkdir(tax_term);
                end
                cd(tax_term)
            end
            lineage5 = unblank(info.organism.lineage{5});
            sfid = fopen(sprintf('%s_disorder_characteristics.csv',lineage5),'a');
            fprintf(sfid,'%s,%5.3f,%6.3f,%6.3f,%i,%i,%s,%s\n',...
                uniprot_ids{t},fraction,ffuzzy,fresidual,psize,nd,go_ids,taxonomy);            fclose(sfid);
            cd(go_home);
        else
            fprintf(ofid,'%s,%5.3f,%6.3f,%6.3f,%i,%i,%s,%s\n',...
                uniprot_ids{t},fraction,ffuzzy,fresidual,psize,nd,go_ids,taxonomy);
        end
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
if ~options.lineage
    fclose(ofid);
end

characteristics.fIDR = all;
characteristics.fuzziness = fuzziness;
characteristics.residual_structure = residual_structure;
characteristics.sequence_lengths = all_sizes;
characteristics.nIFR = nIFR;
characteristics.proteins = proteins;
characteristics.oversized = oversized;

function term = unblank(term)

for k = 1:length(term)
    if term(k) == ' '
        term(k) = '_';
    end
end