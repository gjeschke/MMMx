function [entity,exceptions] = get_AF(UniProtID,options)
%
% GET_AF Load predicted structure from AlphaFold Protein Structure Database
%
%   [entity,exceptions] = GET_AF(UniProtID)
%   Returns an (entity) structure in MMMx:atomic representation
%
% INPUT
% UniProtID     UniProt identifier, for instance 'Q07955' for SRSF1
%               empty output and exception, when a prediction is not 
%               available for this UniProt entry
% options       .keep_file  if true, the downloaded PDB file is kept,
%                           defaults to false (file is deleted)
%
% OUTPUT
% entity       entity structure in MMMx:atomic representation 
%              .name            first four characters of the UniProt
%                               identifier
%              in addition to PDB entities, AlphaFold entities contain
%              .origin          'AlphaFold v%i', where %i is the version
%                               number
%              .uniprot         UniProt identifier
%              .uniprotname     protein name in UniProt
%              .sequence        amino acid sequence
%              .pae             AlphaFold predicted aligned error matrix
% exceptions   error message if something went wrong, entity is empty for
%              errors but not for warnings, for warnings, only the last one
%              is reported, cell containing an empty array, if no exception
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

% initialize empty outputs
exceptions{1} = [];
entity = [];

if ~exist('options','var') || isempty(options)
    options.keep_file = false;
end

query = sprintf('https://alphafold.ebi.ac.uk/api/prediction/%s',UniProtID);
try
    AF_info = webread(query);
catch exception
    exceptions{1} = exception;
    return;
end

fname = sprintf('AF_%s_v%i.pdb',UniProtID,AF_info.latestVersion);
try
    websave(fname,AF_info.pdbUrl);
catch exception
    exceptions{1} = exception;
    return;
end
entity = get_pdb(fname);
entity.name = AF_info.uniprotAccession(1:4);
entity.uniprot = AF_info.uniprotAccession;
entity.uniprotname = AF_info.uniprotId;
entity.sequence = AF_info.uniprotSequence;
entity.origin = sprintf('AlphaFold v%i',AF_info.latestVersion);

wroptions = weboptions('ContentType','json');
pae = webread(AF_info.paeDocUrl,wroptions);
entity.pae = pae.predicted_aligned_error;

if ~options.keep_file
    delete(fname);
end

