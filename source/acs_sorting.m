function [D,indices,passes,h] = acs_sorting(D)
%
% ACS_SORTING   Sorts conformers by distance in abstract conformer space
%
%   [D,indices] = acs_sorting(D,Rg)
%
% INPUT
% D             dRMSD matrix
%
% OUTPUT
% D             sorted dRMSD matrix
% indices       sorted indices for the original dRMSD matrix
% passes        number of passes
% h             figure handle for a plot, the figure is generated only if
%               the handle is requested
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2025: Gunnar Jeschke

[C,~] = size(D); % number of conformers

weights = zeros(C);
for c1 = 1:C-1
    for c2 = c1+1:C
        weights(c1,c2) = abs(c1-c2);
        weights(c2,c1) = weights(c1,c2);
    end
end
weights = weights.^2;

diff = 1/sum(sum(weights.*D));

tic
indices = 1:C;
maxpasses = 150;
passes = 0;
all_diff = zeros(1,maxpasses*C*(C-1)/2);
poi = 0;
changed = true;
while passes < maxpasses && changed
    passes = passes +1;
    changed = false;
    for c1 = 1:C-1
        for c2 = c1+1:C
            poi = poi + 1;
            indices2 = indices;
            indices2(c1) = indices(c2);
            indices2(c2) = indices(c1);
            D2 = D(indices2,indices2);
            dtest = 1/sum(sum(weights.*D2));
            all_diff(poi) = dtest;
            if dtest < diff
                diff = dtest;
                indices = indices2;
                changed = true;
            end
        end
    end
end
% toc,
% fprintf(1,'Sorting took %i passes\n',pass);

D = D(indices,indices);

if nargout >= 4
    h = plot_pair_drmsd(D);
end

% Diagnostics, display new D

function h = plot_pair_drmsd(D)

h = figure; hold on

image(D,'CDataMapping','scaled');
curr_axis = gca;
curr_axis.YDir = 'normal';
colorbar;
axis tight
xlabel('Conformer number');
ylabel('Conformer number');
title('Sorted dRMSD');
axis equal
