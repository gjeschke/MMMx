function name = mk_AF3_input(UniProtID,nseeds,options)
%
% MK_AF3_INPUT Make AlphaFold3 .json input files from UniProt identifier
%
%   mk_AF3_input(UniProtID,nseeds,options)
%     nseeds input files are written to a subdirectory named by options.name 
%     or UniProtID
%
% INPUT
% UniProtID     UniProt identifier, for instance 'Q07955' for SRSF1
%               no action when the identifier does not exist, if empty, the
%               sequence must be provided in options
% nseeds        number of random seeds, which is also number of input files
% options       additional options, structure with fields
%               .copies     number of protomers in a multimer, default: 1
%               .RNA        an RNA sequence
%               .RNA_copies number of copies of the RNA, defaults to 1
%               .ligands    array of struct with a list of ligands with fields
%                           .ID     identifier, suchs as 'K1'
%                           .CCD    CCD code, such as 'K' for potassium ion
%               .sequence   sequence in single-letter format, used only if
%                           UniProtID is empty
%               .name       
%
% OUTPUT
% name         name of the subdirectory with the files, which is also the
%              basis name of the individual files
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2024: Gunnar Jeschke

if ~exist('options','var') || isempty(options) || ~isfield(options,'copies')
    options.copies = 1;
end

if isfield(options,'RNA') || ~isfield(options,'RNA_copies') || isempty(options.RNA_copies)
    options.RNA_copies = 1;
end

if ~isfield(options,'name') || isempty(options.name)
    if ~isempty(UniProtID)
        options.name = UniProtID;
    else
        options.name = 'AlphaFold3';
    end
end

name = options.name;

go_home = pwd;
mkdir(options.name);
cd(options.name);

if isempty(UniProtID)
        sequence{1} = options.sequence;
        type{1} = 'proteinChain'; % can also be rnaSequence
        copies{1} = options.copies;
        id{1} = '';
else
    if ~iscell(UniProtID)
        seq = get_sequence(UniProtID);
        sequence{1} = seq;
        type{1} = 'proteinChain'; % can also be rnaSequence
        copies{1} = options.copies;
        id{1} = '';
    else
        for p = 1:length(UniProtID)
            seq = get_sequence(UniProtID{p});
            sequence{p} = seq; %#ok<*AGROW> 
            type{p} = 'proteinChain'; % can also be rnaSequence
            copies{p} = options.copies(p);
            id{p} = '';
        end
    end
end

seeds = zeros(nseeds,1);
seeds(1) = round(1e9*rand);
for s = 2:nseeds
    seed = 0;
    unique = min(abs(seeds(1:s-1)-seed));
    while ~unique || seed < 1 || seed > 999999999
        seed = round(1e9*rand);
    end
    seeds(s) = seed;
end

bas = length(sequence);

if isfield(options,'RNA') && ~isempty(options.RNA)
    sequence{bas+1} = options.RNA;
    type{bas+1} = 'rnaSequence';
    copies{bas+1} = options.RNA_copies;
    bas = bas + 1;
end

if isfield(options,'ligands') && ~isempty(options.ligands)
    for k = 1:length(options.ligands)
        sequence{bas+k} = options.ligands(k).CCD; 
        type{bas+k} = 'ligand'; 
        copies{bas+k} = options.ligands(k).copies; 
        id{bas+k} = options.ligands(k).ID;
    end
    bas = bas + length(options.ligands);
end

if isfield(options,'ions') && ~isempty(options.ions)
    for k = 1:length(options.ions)
        sequence{bas+k} = options.ions(k).ID; 
        type{bas+k} = 'ion'; 
        copies{bas+k} = options.ions(k).copies; 
    end
end

for s = 1:nseeds
    fname = sprintf('%s_%i.json',options.name,s);
    fid = fopen(fname,'wt');
    fprintf(fid,'[{\n');
    fprintf(fid,' "name": "%s_%i",\n',options.name,s);
    fprintf(fid,' "modelSeeds": [\n');
    fprintf(fid,'  "%i"\n',seeds(s));
    fprintf(fid,' ],\n');
    fprintf(fid,' "sequences": [\n');
    for seq = 1:length(sequence)
        fprintf(fid,'   {\n');
        fprintf(fid,'    "%s": {\n',type{seq});
        switch type{seq}
            case {'proteinChain','rnaSequence'}
                fprintf(fid,'    "sequence":\n');
                fprintf(fid,'    "%s",\n',sequence{seq});
            case {'ligand'}
                fprintf(fid,'    "ligand": "%s",\n',sequence{seq});
            case {'ion'}
                fprintf(fid,'    "ion": "%s",\n',sequence{seq});
        end
        fprintf(fid,'    "count": %i\n',copies{seq});
        fprintf(fid,'   }\n');
        fprintf(fid,'  }');
        if seq < length(sequence)
            fprintf(fid,',');
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,' ],\n');
    fprintf(fid,'"dialect": "alphafoldserver",\n');
    fprintf(fid,'"version": 1\n');
    fprintf(fid,'}\n');
    fprintf(fid,']\n');
    fclose(fid);
end

cd(go_home);