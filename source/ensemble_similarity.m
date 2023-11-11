function [similarity,width1,width2,dev12,Rg1,Rg2] = ensemble_similarity(entity1,chain1,range1,entity2,chain2,range2)

if ~exist('chain1','var') 
    chain1 = '';
end
if ~exist('range1','var') 
    range1 = [];
end
if ~exist('chain2','var') 
    chain2 = '';
end
if ~exist('range2','var') 
    range2 = [];
end

[pair_drms,~,Rg] = pair_drms_matrix(entity1,chain1,range1,entity2,chain2,range2);
m1 = length(entity1.populations);
all_Rg1 = Rg(1:m1);
Rg1 = sqrt(sum(entity1.populations.*all_Rg1.^2));
all_Rg2 = Rg(m1+1:end);
Rg2 = sqrt(sum(entity2.populations.*all_Rg2.^2));
matches = pair_drms(1:m1,m1+1:end);
selfmatch1 = pair_drms(1:m1,1:m1);
selfmatch2 = pair_drms(m1+1:end,m1+1:end);
var11 = kron(entity1.populations,entity1.populations').*selfmatch1.^2;
var22 = kron(entity2.populations,entity2.populations').*selfmatch2.^2;
var12 = kron(entity1.populations,entity2.populations').*matches.^2;
width1 = sqrt(sum(sum(var11)));
width2 = sqrt(sum(sum(var22)));
dev12 = sqrt(sum(sum(var12)));
similarity = sqrt((width1*width2)/sum(sum(var12)));
