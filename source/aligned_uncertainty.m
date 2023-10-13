function [AU,resnums,entity] = aligned_uncertainty(entity,chain,range)
%
% ALIGNED_UNCERTAINTY   Aligned uncertainty matrix for pairwise relative 
%                       residue (CA) positions of a single chain
%
%   [AU,entity] = aligned_uncertainty(entity,chain,range)
%   Augments the entity, writes a rigid-body PDB file, and the frame of a
%   RigiFlex modelling script
%
%   entity = domain_partitioning(entity,options)
%   Allows for control of partitioning and of output
%
% INPUT
% entity        MMMx:atomic entity
% chain         chain, for which aligned uncertainty is to be computed,
%               defaults to 'A'
% range         double(1,2), optional residue range, default: all residues
%
% OUTPUT
% AU           aligned uncertainty matrix size [R,R], where R is the number
%              of residues in the chain
% resnums      residue numbers
% entity       entity with added fields entity.au.(chain).resnums 
%
% Remarks:  This is similar in spirit to predicted aligned error in
%           AlphaFold, AU(i,j) is the position rmsd of the CA arom of
%           residue j when aligning all conformers at the N, CA, and C
%           atoms of residue i

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

% set defaults for missing input
if ~exist('chain','var') || isempty(chain)
    chain = 'A';
end
if ~exist('range','var') || isempty(range)
    range = [-1e01, 1e10];
end

options.full = true;
[backbones,pop] = get_backbones_ensemble(entity,chain,range,options);

C = length(pop); % number of conformers
resnums = backbones.(chain).mono; % residue numbers
R = length(resnums); % number of residues

AU = zeros(R);
N1 = backbones.(chain).N{1};
CA1 = backbones.(chain).bb{1};
C1 = backbones.(chain).C{1};
all_CAj = cell(1,C);
for r1 = 1:R
    coor1 = [N1(r1,:);CA1(r1,:);C1(r1,:)];
    all_CAj{1} = CA1;
    for c = 2:C
        N = backbones.(chain).N{c};
        CA = backbones.(chain).bb{c};
        CO = backbones.(chain).C{c};
        coor2 = [N(r1,:);CA(r1,:);CO(r1,:)];
        [~, ~, transmat] = superimpose_3points(coor1,coor2);
        all_CAj{c} = affine_coor_set(CA,transmat);
    end
    for r2 = 1:R
        CAj = zeros(C,3);
        for c = 1:C
            coor = all_CAj{c};
            CAj(c,:) = coor(r2,:);
        end
        mean_coor = mean(CAj);
        CAj = CAj - repmat(mean_coor,C,1);
        AU(r1,r2) = sqrt(sum(pop'*CAj.^2));
    end
end


entity.au.(chain).resnums = resnums;
entity.au.(chain).AU = AU;

% AU(AU>31.75) = 31.75; % capping like in AlphaFold

figure; clf; hold on
image(resnums,resnums,AU,'CDataMapping','scaled');
curr_axis = gca;
set(curr_axis,'YDir','normal');
colorbar;
axis tight
xlabel('Residue number');
ylabel('Residue number');
title('Aligned uncertainty');
axis equal