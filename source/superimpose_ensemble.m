function [entity,rmsd,exceptions] = superimpose_ensemble(entity,selected,template)
%
% SUPERIMPOSE_ENSEMBLE Superimposes conformers in an ensemble
%
%   [entity,exceptions] = SUPERIMPOSE_ENSEMBLE(entity)
%   Optimal superposition of all coordinates onto first conformer 
%
%   [entity,exceptions] = SUPERIMPOSE_ENSEMBLE(entity,selected)
%   Optimal superposition at selected chains, residues, or atoms
%
%   [entity,exceptions] = SUPERIMPOSE_ENSEMBLE(entity,selected,template)
%   Optimal superposition at selected chains, residues, or atoms to a
%   template structure
%
% INPUT
% entity    entity that describes the ensemble
% selected  MMMx address of the atoms that are to be superimposed, defaults
%           to all atoms
%           string: the same address is used in entity and template
%           cell string containing two strings: the first string is used in
%           entity and the second string in template
% template  template entity, defaults to none
%
% OUTPUT
% entity       entity after superposition, only coordinates have changed
% exceptions   error message if something went wrong, entity is empty for
%              errors but not for warnings, for warnings, only the last one
%              is reported, cell containing an empty array, if no exception
% rmsd         vector of roor mean square atom deviations for all
%              conformers
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

rmsd = zeros(length(entity.populations),1);

if ~exist('selected','var')
    selected = '(*)';
end

if iscell(selected)
    selection_template = selected{2};
    selected = selected{1};
else
    selection_template = selected;
end

if ~exist('template','var') || isempty(template)
    [tcoor,~,exceptions] = get_coor(entity,['{1}' selected]);
else
    [tcoor,~,exceptions] = get_coor(template,selection_template);
end

if ~isempty(exceptions) 
    if ~isempty(exceptions{1})
        return
    end
end

for c = 1:length(entity.populations)
    [ccoor,~,exceptions] = get_coor(entity,[sprintf('{%i}',c) selected]);
    if ~isempty(exceptions)
        if ~isempty(exceptions{1})
            return
        end
    end
    [rms,~,transmat] = rmsd_superimpose(tcoor,ccoor);
    rmsd(c) = rms;
    [ccoor,indices,exceptions] = get_coor(entity,sprintf('{%i}(*)',c));
    if ~isempty(exceptions)
        if ~isempty(exceptions{1})
            return
        end
    end
    xyz = [ccoor ones(length(indices),1)]*transmat';
    entity.xyz(indices,:) = xyz(:,1:3);
end
