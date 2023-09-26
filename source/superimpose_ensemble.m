function [entity,rmsd,exceptions] = superimpose_ensemble(entity,selected,template,central,options)
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
%           string: the same address is used in entity and template, if
%                   empty, whole structure is selected
%           cell string containing two strings: the first string is used in
%           entity and the second string in template
% template  template entity, defaults to none
% central   flag that requests superposition onto the central structure
%           if a template is provided, superposition is onto the common
%           central structure, defaults to false
%
% OUTPUT
% entity       entity after superposition, only coordinates have changed
% exceptions   error message if something went wrong, entity is empty for
%              errors but not for warnings, for warnings, only the last one
%              is reported, cell containing an empty array, if no exception
% rmsd         vector of root mean square atom deviations for all
%              conformers
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

if ~exist('options','var') || isempty(options)
    options.atoms = 'all';
end

rmsd = zeros(length(entity.populations),1);

if ~exist('selected','var') || isempty(selected)
    selected = '(*)';
end

if iscell(selected)
    selection_template = selected{2};
    selected = selected{1};
else
    selection_template = selected;
end

switch options.atoms
    case 'backbone'
        selection_template = strcat(selection_template,'.N,CA,C,O');
        selected = strcat(selected,'.N,CA,C,O');
    case 'CA'
        selection_template = strcat(selection_template,'.CA');
        selected = strcat(selected,'.CA');
end
if ~exist('template','var') || isempty(template)
    [tcoor,~,~,exceptions] = get_coor(entity,['{1}' selected]);
else
    [tcoor,~,~,exceptions] = get_coor(template,selection_template);
end

if ~exist('central','var') || isempty(central)
    central = false;
end

% superposition to central structure overrides template selection or
% selection in first conformer
if central
    [chain1,range1] = split_chain_range(selected);
    if ~exist('template','var') || isempty(template)
        tcoor = central_structure_coor(entity,chain1,range1);
    else
        tcoor = central_structure_coor(entity,selected,template,selection_template);
    end
end

if ~isempty(exceptions) 
    if ~isempty(exceptions{1})
        return
    end
end

for c = 1:length(entity.populations)
    [ccoor,~,~,exceptions] = get_coor(entity,[sprintf('{%i}',c) selected]);
    % [ccoor,~,~,exceptions] = get_coor(entity,[sprintf('{%i}',c) selected '.N,CA,C,O']);
    if ~isempty(exceptions)
        if ~isempty(exceptions{1})
            return
        end
    end
    [rms,~,transmat] = rmsd_superimpose(tcoor,ccoor);
    rmsd(c) = rms;
    indices = find(entity.index_array(:,4) == c);
    ccoor = entity.xyz(indices,:);
    xyz = [ccoor ones(length(indices),1)]*transmat';
    entity.xyz(indices,:) = xyz(:,1:3);
end

function coor = central_structure_coor(entity1,selected1,entity2,selected2)
% for selecting everything, chain and range identifiers should be empty

C = length(entity1.populations);
if exist('entity2','var') && ~isempty(entity2) && exist('selected2','var')
    C2 = length(entity2.populations);
else
    C2 =0;
end
all_coor = cell(1,C+C2);
for c = 1:C
    % combine coordinate arrays
    coor = get_coor(entity1,[sprintf('{%i}',c) selected1]);
    all_coor{c} = coor;
end

for c = C+1:C+C2
    % combine coordinate arrays
    coor = get_coor(entity2,[sprintf('{%i}',c-C) selected2]);
    all_coor{c} = coor;
end
C = C + C2;

pair_rmsd = zeros(C);
for c1 = 1:C-1
    coor1 = all_coor{c1};
    for c2 = c1+1:C
        coor2 = all_coor{c2};
        rmsd = rmsd_superimpose(coor1,coor2);
        pair_rmsd(c1,c2) = rmsd;
        pair_rmsd(c2,c1) = rmsd;
    end
end

sum_msq = sum(pair_rmsd.^2);
[~,conf] = min(sum_msq);
coor = all_coor{conf};


% function coor = central_structure_coor_bb(entity1,chain1,range1,entity2,chain2,range2)
% % for selecting everything, chain and range identifiers should be empty
% 
% backbones = get_backbones_ensemble(entity1,chain1,range1);
% chains = fieldnames(backbones);
% N = 0;
% for kc = 1:length(chains)
%     bb0 = backbones.(chains{kc}).bb{1};
%     [N0,~] = size(bb0);
%     N = N + N0;
% end
% backbones.all.chains = blanks(N);
% backbones.all.mono = zeros(1,N);
% backbones.all.slc = blanks(N);
% poi = 0;
% for kc = 1:length(chains)
%     for km = 1:length(backbones.(chains{kc}).mono)
%         backbones.all.chains(poi+km) = chains{kc};
%         backbones.all.slc(poi+km) = backbones.(chains{kc}).slc(km);
%         backbones.all.mono(poi+km) = backbones.(chains{kc}).mono(km);
%     end
%     poi = poi + length(backbones.(chains{kc}).mono);
% end
% C = length(backbones.(chains{1}).bb);
% for c = 1:C
%     % combine coordinate arrays
%     bb = zeros(N,3);
%     poi = 0;
%     for k = 1:length(chains)
%         bb0 = backbones.(chains{k}).bb{c};
%         [N0,~] = size(bb0);
%         bb(poi+1:poi+N0,:) = bb0;
%         poi = poi + N0;
%     end
%     % assign complete backbone
%     backbones.all.bb{c} = bb;
% end
% 
% if exist('entity2','var') && ~isempty(entity2) && exist('chain2','var') && exist('range2','var')
%     backbones2 = get_backbones_ensemble(entity2,chain2,range2);
%     chains2 = fieldnames(backbones2);
%     N = 0;
%     for kc = 1:length(chains2)
%         bb0 = backbones2.(chains2{kc}).bb{1};
%         [N0,~] = size(bb0);
%         N = N + N0;
%     end
%     C2 = length(backbones.(chains{1}).bb);
%     for c = C+1:C+C2
%         % combine coordinate arrays
%         bb = zeros(N,3);
%         poi = 0;
%         for k = 1:length(chains)
%             bb0 = backbones.(chains{k}).bb{c-C};
%             [N0,~] = size(bb0);
%             bb(poi+1:poi+N0,:) = bb0;
%             poi = poi + N0;
%         end
%         % assign complete backbone
%         backbones.all.bb{c} = bb;
%     end
%     C = C + C2;
% end
% 
% 
% pair_rmsd = zeros(C);
% for c1 = 1:C-1
%     coor1 = backbones.all.bb{c1};
%     for c2 = c1+1:C
%         coor2 = backbones.all.bb{c2};
%         rmsd = rmsd_superimpose(coor1,coor2);
%         pair_rmsd(c1,c2) = rmsd;
%         pair_rmsd(c2,c1) = rmsd;
%     end
% end
% 
% sum_msq = sum(pair_rmsd.^2);
% [~,conf] = min(sum_msq);
% coor = backbones.all.bb{conf};
% 
function [chain,range] = split_chain_range(address)

range = [];

if contains(address,'*')
    chain = '';
    return;
end

poia = strfind(address,'(');
poie = strfind(address,')');
if isempty(poie)
    chain = '';
else
    chain = address(poia+1:poie-1);
    address = address(poie+1:end);
end
if strcmp(chain,'*')
    range = [1,1e6];
    return
end
residues = split(address,'-');
if ~isempty(residues{1})
    range(1) = str2double(residues{1});
end
if ~isempty(residues{2})
    range(2) = str2double(residues{2});
end
if length(range) == 1
    range(2) = range(1);
end