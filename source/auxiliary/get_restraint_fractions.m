function fractions = get_restraint_fractions(tindex,intervals,digitbases)
% fractions = get_restraint_fractions(tindex,digitbases)
%
% Computes the fractions of the lower-bound to upper-bound range for all
% restraints that correspond to a trial index
% (subroutine for the exhaustive sampling Rigi module of RigiFlex)
%
% tindex        trial index
% intervals     number of intervals per restraint, vector of length
%               9*r*(r-1)/2 for r rigid bodies
% digitbases    number of trials for earlier restraints, vector
%               [1,9*r*(r-1)/2] 
%
% fractions     range fractions, vector [1,9*r*(r-1)/2]
%
% G. Jeschke, 3.10.2017

m = length(intervals);
fractions = zeros(1,m);
tindex = tindex-1;
for k = m:-1:1
    digit = floor(tindex/digitbases(k));
    tindex = tindex - digit*digitbases(k);
    fractions(k) = (digit + 1)/(intervals(k)+1);
end