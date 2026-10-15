function flexibility = flexibility_parameters_protein(UniProtId,GPR_disorder,GPR_RMSF,GPR_SSP)
% FLEXIBILITY_PARAMETERS_proteome(UPID_file,options)
% Determines some disorder characteristics for a list of UniProt
% identifiers and writes a file with the characteristics for individual
% proteins as well as files for all gene ontology identifiers
%
% Input:
% UniProtId     UniProt identifier of the protein
% GPR_disorder  Gaussian process regeression model for disorder prediction,
%               optional, is loaded if missing
% GPR_RMSF      GPR model for predicting root mean square fluctuation,
%               optional, is loaded if missing
% GPR_SSP       GPR model for predicting secondary structure propensity
%               optional, is loaded if missing
%
% G. Jeschke, 2026

if ~exist('GPR_disorder','var')
    GPR_disorder = load('predict_disorder_RFM_from_AF.mat');
end
if ~exist('GPR_RMSF','var')
    GPR_RMSF = load('predict_RMSF_RFM_from_AF.mat');
end
if ~exist('GPR_SSP','var')
    GPR_SSP = load('predict_SSP_RFM_from_AF.mat');
end

AF_options.structure = false;
AF_options.pLDDT = true;
AF_options.pae = true;

entity = get_AF(UniProtId,AF_options);
if isempty(entity)
    fprintf(2,'No AlphaFold data for UniProtId %s. Aborting.\n',UniProtId);
    flexibility = [];
    return
end
flexibility = get_flexibility(entity,GPR_disorder,GPR_RMSF,GPR_SSP);

