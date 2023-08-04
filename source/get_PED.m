function [entity,exceptions] = get_PED(pedID,ensembleID,options)
%
% GET_PED Load ensemble structure from Protein Ensemble Database (PED)
%
%   [entity,exceptions] = GET_PED(pedID,ensembleID,options)
%   Returns an (entity) structure in MMMx:atomic representation
%
% INPUT
% pedID         PED identifier, for instance 'ped00057' for PARP-1/DNA
%               empty output and exception, when this PED ID does not
%               exist
% ensembleID    ensemble identifier, for instance 'e001', empty output and
%               exception, when this ensemble ID does not exist
% options       .keep_file  if true, the downloaded and extracted PDB file
%                           is kept, defaults to false (file is deleted)
%
% OUTPUT
% entity       entity structure in MMMx:atomic representation 
%              .name            E, followed by last three digits of PED ID
%              in addition to PDB entities, PED entities contain
%              .title           title of the PED entry
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

if ~options.keep_file
    delete(fname);
end

