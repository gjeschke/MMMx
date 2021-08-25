function alignment = sequence_alignment(seqs)
% alignment = SEQUENCE_ALIGNMENT(seqs)
%
% aligns peptide sequences by calling MUSCLE by Robert C. Edgar
% citation: Edgar, R.C. Nucleic Acids Res 32(5), 1792-97
%
% INPUT:
%
% seqs          cell array of peptide sequences (at least two) 
%
% OUTPUT:
%
% alignment     array[s,N] of numbers of corresponding residues, where
%               numbers start with 1 for each of the s sequences, NaN
%               N corresponds to the total length of the alignment 
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% Initialuize empty output
alignment = [];

for k = 1:length(seqs) % replace unknown residue identifier to conform to common usage
    seq = seqs{k};
    for kk = 1:length(seq)
        if char(seq(kk)) == '?'
            seq(kk) = 'X';
        end
    end
    seqs{k} = seq;
end

if length(seqs) < 2 % return with empty output if only one or none sequence was provided
    return
end

alignment(length(seqs)).db = 'MMMx';
alignment(length(seqs)).name = '';
alignment(length(seqs)).sequence = '';

for k = 1:length(seqs)
    alignment(k).db = 'MMMx';
    alignment(k).name= sprintf('seq%i',k);
    alignment(k).sequence = seqs{k};
end

dospath = which('muscle.exe');
if isempty(dospath) % return with empty alignment if MUSCLE is missing
    return
end

infile = fullfile(pwd,'doalign.fa');
wr_multiple_fasta(infile,alignment);
outfile = fullfile(pwd,'muscle_alignment.aln');
comd = [dospath ' -in ' infile ' -clwout ' outfile];
if ispc
    [s,~] = dos(comd);
elseif isunix
    [s,~] = unix(comd);
end
if s ~= 0 % return with empty output, if MUSCLE failed
    return
end

alignment_records = rd_multiple_clustal(outfile);

if ~isempty(alignment_records)
    N = length(alignment_records(1).sequence);
    s = length(seqs);
    alignment = NaN(s,N);
    counters = zeros(s,1);
    for kr = 1:N
        for ks = 1:s
            if alignment_records(ks).sequence(kr) ~= '-'
                counters(ks) = counters(ks) + 1;
                alignment(ks,kr) = counters(ks);
            end
        end
    end
end

delete(outfile);
delete(infile);


function wr_multiple_fasta(fname,alignment)
% function wr_multiple_fasta(fname,alignment)
%
% generates a multiple sequence file in FASTA format when given a multiple 
% sequence record or alignment and, optionally, a list of
% sequence numbers to be included, this can be used, among else, for
% creating input for the sequence alignment program MUSCLE by R. C. Edgar
%
% alignments can be read with 
% get_multiple_fasta.m
% get_multiple_pir.m
% get_multiple_clustal.m
%
% fname     full file name, including extension
% alignment array of alignment struct, as created by the read routines mentioned above
%           as a minimum, each item must have fields
%           .db
%           .sequence
%           .name           
%
% Remarks:
% - assumed format for the so-called FASTA defline is based on the 
%   "Sequence identifiers" section of
%   http://en.wikipedia.org/wiki/FASTA_format (as of 21.03.2011)
%
%   generally, the FASTA defline is ill-defined, but major databases use
%   well-defined subformats, for pdb we use the NCBI Blast flavor of the
%   defline, rather than the very non-standard format provided by the PDB
%   homepage
%   sister routine get_multiple_fasta.m accepts both PDB defline formats
%
% G. Jeschke, 2011

if isempty(alignment)
    return
end

fid = fopen(fname,'wt');
if fid == -1
    return;
end

for cseq = 1:length(alignment)
    alg = alignment(cseq);
    % defline preamble, which is independent of database (replacing the
    % messy PDB FASTA format by a sensible one, as NCBI does)
    fprintf(fid,'>%s|',alg.db);
    % rest of defline depends on data base
    switch alg.db
        case 'gi'
            fprintf(fid,'%s|%s|%s| %s\n',alg.gi_number,alg.type,alg.accession,alg.name);
        case 'pir'
            fprintf(fid,'|%s\n',alg.name);
        case 'prf'
            fprintf(fid,'|%s\n',alg.name);
        case 'sp'
            fprintf(fid,'%s|%s\n',alg.accession,alg.name);
        case 'pdb'
            fprintf(fid,'%s|%s\n',alg.name,alg.chain);
        case 'pat'
            fprintf(fid,'%s|%s\n',alg.country,alg.name);
        case 'bbs'
            fprintf(fid,'%s\n',alg.name);
        case 'gnl'
            fprintf(fid,'%s|%s\n',alg.type,alg.name);
        case 'ref'
            fprintf(fid,'%s|%s\n',alg.accession,alg.name);
        case 'tr'
            fprintf(fid,'%s|%s\n',alg.accession,alg.name);
        case 'lcl'
            fprintf(fid,'%s\n',alg.name);
        otherwise
            fprintf(fid,'%s\n',alg.name);
    end
    % write sequence, up to 70 residues per line
    tsequence = alg.sequence;
    while length(tsequence) > 70
        fprintf(fid,'%s\n',tsequence(1:70));
        tsequence = tsequence(71:end);
    end
    fprintf(fid,'%s\n\n',tsequence);
end

fclose(fid);

function alignment = rd_multiple_clustal(inname)
%
% function alignment = rd_multiple_fasta(inname),
%
% Reads multiple peptide (protein) sequences in CLUSTAL format
%
% inname    filename INCLUDING extension
% alignment alignment array of struct (one for each sequence) with fields 
%           .db         data base type
%           .sequence   string that codes the sequence in single-letter format
%           .name       identifier
%           further fields depend on data base type (see Reamrks)
%
% Remarks:
% - assumed input format for the so-called FASTA defline is based on the 
%   "Sequence identifiers" section of
%   http://en.wikipedia.org/wiki/FASTA_format (as of 21.03.2011)
%
%   generally, the FASTA defline is ill-defined, but major databases use
%   well-defined subformats
%
%
% (c) G. Jeschke, 2011

alignment(100).db = 'Clustal';
alignment(100).name = '';
alignment(100).sequence = '';

rfile = fopen(inname,'r');

taglist = '\';
nl=0;
idmax =0;
while 1
    tline = fgetl(rfile); nl = nl+1;
    if ~ischar(tline), break, end
    if nl == 1
        if ~strcmpi(tline(1:7),'CLUSTAL') && ~strcmpi(tline(1:6),'MUSCLE')
            break
        end
        tline=fgetl(rfile); nl=nl+1;
        if ~ischar(tline), break, end
    end
    if ~isempty(strtrim(tline))
        nonsense = textscan(tline,'%s');
        if char(tline(1)) ~= ' ' && ~isempty(nonsense) && length(nonsense{1}) >= 2 % line with sequence information
            tags = nonsense{1};
            id = tag2id(char(tags(1)),taglist,[],'\'); % find sequence number
            if isempty(id)
                taglist = [taglist char(tags(1)) '\']; %#ok<AGROW> % add tag, if sequence was still unknown
                id = tag2id(char(tags(1)),taglist,[],'\');
                alignment(id).db = 'Clustal';
                alignment(id).name = char(tags(1));
                alignment(id).sequence = '';
            else
                if id > idmax
                    idmax = id;
                end
            end
            alignment(id).sequence = [alignment(id).sequence char(tags(2))];
        end
    end
end
fclose(rfile);

alignment = alignment(1:idmax);