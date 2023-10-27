function [entity,exceptions] = get_PED(pedID,ensembleID,options)
%
% GET_PED Load ensemble structure from Protein Ensemble Database (PED)
%
%   [entity,exceptions] = GET_PED(pedID,ensembleID,options)
%   Returns an (entity) structure in MMMx:atomic representation
%
% INPUT
% pedID         PED identifier, for instance 'PED00057' for PARP-1/DNA
%               empty output and exception, when this PED ID does not
%               exist
% ensembleID    (optional) ensemble identifier, for instance 'e001', 
%               empty output and exception, when this ensemble ID does not 
%               exist, if input is missing, all ensembles are downloaded
% options       .keep_file  if true, the downloaded and extracted PDB file
%                           is kept, defaults to false (file is deleted)
%
% OUTPUT
% entity       entity structure in MMMx:atomic representation 
%              if input 'ensembleID' is missing, cell of entities with all
%              ensembles in this PED entry
%              .name            E, followed by last three digits of PED ID
%              in addition to PDB entities, PED entities contain
%              .title           title of the PED entry
%              .ensembleID             ensemble IF
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

query = sprintf('https://deposition.proteinensemble.org/api/v1/entries/%s',pedID);
try
    PED_info = webread(query);
catch exception
    exceptions{1} = exception;
    return;
end

if ~exist('ensembleID','var')
    eIDs = cell(1,length(PED_info.ensembles));
    for ens = 1:length(PED_info.ensembles)
        eIDs{ens} = PED_info.ensembles(ens).ensemble_id;
    end
else
    eIDs{1} = ensembleID;
end

entities = cell(1,length(eIDs));

for ens = 1:length(eIDs)
    
    ensembleID = eIDs{ens};
    query = sprintf('https://deposition.proteinensemble.org/api/v1/entries/%s/ensembles/%s/ensemble-pdb',pedID,ensembleID);
    fname = sprintf('%s_%s.tar.gz',pedID,ensembleID);
    tarname = sprintf('%s_%s.tar',pedID,ensembleID);
    try
        websave(fname,query);
    catch exception
        exceptions{1} = exception;
        return;
    end
    gunzip(fname);
    delete(fname);
    filenames = untar(tarname);
    delete(tarname);
    fname = filenames{1};
    entity = get_pdb(fname);
    entity.title = PED_info.description.title;
    entity.name = sprintf('E%s',pedID(end-2:end));
    entity.ensembleID = ensembleID;
    entities{ens} = entity;
    if ~options.keep_file
        delete(fname);
    end
    
    
    query = sprintf('https://deposition.proteinensemble.org/api/v1/entries/%s/ensembles/%s/weights',pedID,ensembleID);
    fname = sprintf('%s_%s_weights.dat',pedID,ensembleID);
    try
        websave(fname,query);
    catch exception
        exceptions{1} = exception;
        return;
    end
    
    weights = load(fname);
    entity.populations = weights(:,2);
    
    if ~options.keep_file
        delete(fname);
    end
    
end

% if more than one ensemble was read, return cell array of all ensembles
if length(eIDs) > 1
    entity = entities;
end