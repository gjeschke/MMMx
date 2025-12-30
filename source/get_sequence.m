function [sequence,header] = get_sequence(UniProtID)
%
% GET_SEQUENCE Get protein sequence from UniProt identifier
%
%   [sequence,header] = get_sequence(UniProtID)
%   Returns the sequence in single-letter format and the header
%
% INPUT
% UniProtID     UniProt identifier, for instance 'Q07955' for SRSF1
%               empty output and exception, when the identifier does not
%               exist
%
% OUTPUT
% sequence     sequence of the canonical form in single-letter format 
% header       FASTA header
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2024: Gunnar Jeschke


query = sprintf('http://www.ebi.ac.uk/Tools/dbfetch/dbfetch/uniprot/%s/fasta',UniProtID);
try
    websave('fasta.seq',query);
    [header,sequence] = get_fasta('fasta.seq');
    delete('fasta.seq');
catch
    sequence = '';
    header = '';
    return;
end
