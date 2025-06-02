function [elDDT,all_lDDT,lEAU,lCI95] = lDDT(entity,rentity,options)
%
% LDDT Compute local distance difference test (lDDT) for an ensemble
%               structure with respect to a reference structure
%
%   elDDT = LDDT(entity,rentity,options)
%   Computes lDDT according to V. Mariani et al., Bioinformatics, 2013, 29,
%   2722-2729 (doi:10.1093/bioinformatics/btt473) for an ensemble with
%   respect to a single reference structure, control parameters can be
%   adjusted via options, 
%   computation of local ensemble aligned uncertainty (EAU) by providing
%   an EAU matrix
%
% INPUT
% entity        entity that contains the ensemble structure
% rentity       entity that contains the reference structure
% options       options, struct
%               .R0             inclusion radius, defaults to 15 Å, connot
%                               be smaller than 4
%               .thresholds     set of tolerance thresholds, defaults to
%                               [0.5, 1, 2, 4] Å
%               .chain          chain identifier, defaults to 'A'
%               .eau            optional EAU matrix, dimension must fit
%                               number of residues
%
% OUTPUT
% elDDT         ensemble lDDT, weighted average over all conformers, range
%               is 0 to 1
% all_lDDT      array (C,R) of lDDT for C conformers with R residues each
% lEAU          local EAU per residue, computed only if options.eau is
%               provided, otherwise empty
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2025: Gunnar Jeschke

% set default values
R0 = 15;
thresholds = [0.5,1,2,4];
chain = 'A';

% check for input options
if exist('options','var')
    if isfield(options,'R0') && ~isempty(options.R0)
        R0 = options.R0;
    end
    if isfield(options,'thresholds') && ~isempty(options.thresholds)
        thresholds = options.thresholds;
    end
    if isfield(options,'chain') && ~isempty(options.chain)
        chain = options.chain;
    end
    if isfield(options,'eau') && ~isempty(options.eau)
        mk_eau = true;
    else
        mk_eau = false;
        lEAU = [];
    end
    if isfield(options,'CI95') && ~isempty(options.CI95)
        mk_CI95 = true;
    else
        mk_CI95 = false;
        lCI95 = [];
    end
end

if R0 < 4 % include at least direct neighbors to avoid division by zero
    R0 = 4;
end

% Determine CA-CA distance matrix of reference structure
resfields = fieldnames(rentity.(chain));
coor = zeros(length(resfields),3);
residues = 0;
for k = 1:length(resfields)
    resstr = resfields{k};
    if resstr(1) == 'R' % this is a residue field
        residues = residues + 1;
        idx = rentity.(chain).(resstr).CA.tab_indices;
        if length(idx) > 1 % if the reference structure is an ensemble, the first conformer acts as the reference
            idx = idx(1);
        end
        coor(residues,:) = rentity.xyz(idx,:);
    end
end
coor = coor(1:residues,:);
ref_D = coor2dmat(coor);

all_lDDT = zeros(length(entity.populations),residues);
if mk_eau
    lEAU = zeros(1,length(residues));
end
if mk_CI95
    lCI95 = zeros(1,length(residues));
end

% determine the lDDT vectors for all conformers
for c = 1:length(entity.populations)
    resfields = fieldnames(entity.(chain));
    coor = zeros(length(resfields),3);
    cresidues = 0;
    for k = 1:length(resfields)
        resstr = resfields{k};
        if resstr(1) == 'R' % this is a residue field
            cresidues = cresidues + 1;
            idx = entity.(chain).(resstr).CA.tab_indices;
            if length(idx) > 1 % if the reference structure is an ensemble, the first conformer acts as the reference
                idx = idx(c);
            end
            coor(cresidues,:) = entity.xyz(idx,:);
        end
    end
    coor = coor(1:cresidues,:);
    D = coor2dmat(coor);
    clDDT = zeros(1,residues);
    for res = 1:residues
        v = ref_D(res,:); % this row of the reference distance matrix
        % find the indices of all residues within the inclusion radius
        indices = find(v > 0.1 & v <= R0);
        n = length(indices); % number of residues in inclusion radius
        % make distance difference vector, absolute numbers
        diff = abs(v(indices) - D(res,indices));
        % find fractions for all threshold values
        fractions = zeros(1,length(thresholds));
        for t = 1:length(thresholds)
            fractions(t) = sum(diff <= thresholds(t))/n;
        end
        clDDT(res) = mean(fractions);
        if mk_eau
            lEAU(res) = mean(options.eau(res,indices)); %#ok<AGROW> 
        end
        if mk_CI95
            lCI95(res) = mean(options.CI95(res,indices)); %#ok<AGROW> 
        end
    end 
    all_lDDT(c,:) = clDDT;
end
elDDT = entity.populations'*all_lDDT;

