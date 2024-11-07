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
    seq = options.sequence;
else
    seq = get_sequence(UniProtID);
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
sequence{1} = seq;
type{1} = 'proteinChain'; % can also be rnaSequence
copies{1} = options.copies;

if isfield(options,'RNA') && ~isempty(options.RNA)
    sequence{2} = options.RNA;
    type{2} = 'rnaSequence';
    copies{2} = options.RNA_copies;
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
        fprintf(fid,'    "sequence":\n');
        fprintf(fid,'    "%s",\n',sequence{seq});
        fprintf(fid,'    "count": %i\n',copies{seq});
        fprintf(fid,'   }\n');
        fprintf(fid,'  }');
        if seq < length(sequence)
            fprintf(fid,',');
        end
        fprintf(fid,'\n');
    end
    fprintf(fid,' ]\n');
    fprintf(fid,'}\n');
    fprintf(fid,']\n');
    fclose(fid);
end

cd(go_home);