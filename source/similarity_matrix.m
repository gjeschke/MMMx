function S = similarity_matrix(D,populations)
% SIMILARITY_MATRIX   Given a distance root mean square deviation matrix
%                     and a vector of conformr numbers in subensembles,
%                     the function computes the similarities between
%                     subensembles
%
%    S = similarity_matrix(dRMSD,C)
% 
% INPUT
% D             distance root mean square deviation matrix for all
%               conformers, must be sorted by subensembles
% populations   cell(N,1) of population vectors of the subensembles
%
% OUTPUT
% S             symmetric square matrix of similarities for all pairs of
%               subensembles, diagonal elements are all ones
%
% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% (c) G. Jeschke, 2025

C = zeros(length(populations),1);
for ens = 1:length(populations)
    C(ens) = length(populations{ens});
end

S = ones(length(C));

indices = zeros(length(C),2);
pointer = 0;
for c = 1:length(C)
    indices(c,1) = pointer + 1;
    indices(c,2) = pointer + C(c);
    pointer = pointer + C(c);
end

for e1 = 1:length(C)-1
    selfmatch1 = D(indices(e1,1):indices(e1,2),indices(e1,1):indices(e1,2));
    pop1 = populations{e1};
    var11 = kron(pop1,pop1').*selfmatch1.^2;
    for e2 = e1+1:length(C)
        selfmatch2 = D(indices(e2,1):indices(e2,2),indices(e2,1):indices(e2,2));
        matches = D(indices(e1,1):indices(e1,2),indices(e2,1):indices(e2,2));
        pop2 = populations{e2};
        var22 = kron(pop2,pop2').*selfmatch2.^2;
        var12 = kron(pop1,pop2').*matches.^2;
        width1 = sqrt(sum(sum(var11)));
        width2 = sqrt(sum(sum(var22)));
        S(e1,e2) = sqrt((width1*width2)/sum(sum(var12)));
        S(e2,e1) = S(e1,e2);
    end
end

