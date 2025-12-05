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
%               .structure  if true, PDB is read for 3D coordinates,
%                           defaults to true
%               .pLDDT      if true, pLDDT is read, defaults to true
%               .pae        if true, pae is read, defaults to true
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
%              .pLDDT           predicted local-distance difference test,
%                               a residue-wise confidence predictor
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

if ~exist('options','var') || isempty(options) || ~isfield(options,'keep_file')
    options.keep_file = false;
end

if ~isfield(options,'structure')
    options.structure = true;
end

if ~isfield(options,'pLDDT')
    options.pLDDT = true;
end

if ~isfield(options,'pae')
    options.pae = true;
end


query = sprintf('https://alphafold.ebi.ac.uk/api/prediction/%s',UniProtID);
try
    AF_info = webread(query);
catch exception
    exceptions{1} = exception;
    return;
end

if iscell(AF_info)
    AF_info = AF_info{1};
else
    AF_info = AF_info(1);
end

fname = sprintf('%s_v%i.pdb',AF_info.entryId,AF_info.latestVersion);
if options.structure
    try
        websave(fname,AF_info.pdbUrl);
    catch exception
        exceptions{1} = exception;
        return;
    end
    entity = get_pdb(fname);
end
entity.name = AF_info.uniprotAccession(1:4);
entity.uniprot = AF_info.uniprotAccession;
entity.uniprotname = AF_info.uniprotId;
entity.sequence = AF_info.uniprotSequence;
entity.origin = sprintf('AlphaFold v%i',AF_info.latestVersion);
entity.organism = AF_info.organismScientificName;
entity.AF_info = AF_info;

wroptions = weboptions('ContentType','json');
if options.pae
    pae = webread(AF_info.paeDocUrl,wroptions);
    entity.pae = pae.predicted_aligned_error;
end
if options.pLDDT
    pLDDT = webread(AF_info.plddtDocUrl,wroptions);
    entity.pLDDT = pLDDT.confidenceScore;
end

if ~options.keep_file && exist(fname,'file')
    delete(fname);
end

% compile pLDDT information

% pLDDT = zeros(1,length(entity.sequence));
% chain = 'A';
% residues = fieldnames(entity.(chain));
% rpoi = 0;
% 
% for kr = 1:length(residues) % expand over all residues
%     residue = residues{kr};
%     if strcmp(residue(1),'R') % these are residue fields
%         rpoi = rpoi + 1;
%         pLDDT(rpoi) = entity.(chain).(residue).CA.bfactor;
%     end
% end
% 
% entity.pLDDT = pLDDT;