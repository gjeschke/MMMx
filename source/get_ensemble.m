function [entity,exceptions] = get_ensemble(fname,name)
%
% GET_ENSEMBLE Generate entity from an MMMx ensemble list
%
%   [entity,exceptions] = GET_ENSEMBLE(fname)
%   Returns an (entity) structure in MMMx:atomic representation
%
%   [entity,exceptions] = GET_ENSEMBLE(fname,name)
%   Returns an (entity) structure in MMMx:atomic representation and names
%   it
%
% INPUT
% fname     file name of the ensemble list, extension .ens is assumed if
%           none is given
% name      optional name for the entity, defaults to PDB identifier of
%           first conformer if header line exists; MMMx otherwise
%
% OUTPUT
% entity       entity structure in MMMx:atomic representation
% exceptions   error message if something went wrong, entity is empty for
%              errors but not for warnings, for warnings, only the last one
%              is reported, cell containing an empty array, if no exception
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty outputs
entity = [];

if ~exist('name','var') || isempty(name)
    name = 'MMMx';
end

[all_files,pop,exceptions] = rd_ensemble_definition(fname);

if isempty(pop)
    return
end

[entity,exceptions] = get_pdb(all_files(1).name);

if ~isempty(exceptions) && ~isempty(exceptions{1})
    return
end

for c = 2:length(pop)
    [entity,exceptions] = get_pdb(all_files(c).name,[],entity);
    if ~isempty(exceptions) && ~isempty(exceptions{1})
        return
    end
end
entity.populations = pop;

entity.name = name;