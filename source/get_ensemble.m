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
% REMARKS
% populations are normalized to unity sum
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty outputs
entity = [];
my_dir = pwd;

[all_files,pop,exceptions] = rd_ensemble_definition(fname);

[fpath,~] = fileparts(fname);
if ~isempty(fpath)
    cd(fpath);
end

if isempty(pop)
    cd(my_dir);
    return
end

filenames = cell(1,length(all_files));
for conf = 1:length(all_files)
    filenames{conf} = all_files(conf).name;
end
[entity,exceptions] = entity_from_filelist(filenames);

if ~isempty(exceptions) && ~isempty(exceptions{1})
    cd(my_dir);
    return
end

entity.populations = pop/sum(pop);

if exist('name','var')
    entity.name = name;
elseif isempty(entity.name)
    entity.name = 'MMMx';
end

cd(my_dir);