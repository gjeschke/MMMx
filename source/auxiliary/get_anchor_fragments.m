function codes = get_anchor_fragments(ptac,direction,shortfrag)
% codes = get_anchor_fragments(ptac,direction,shortfrag)
%
% provide indices (codes) of fragments from a nucleotide library that fit
% within a threshold to existing pseudo-torsion defining atom coordinates
%
% ptac      pseudo-torsion defining atom coordinates
% direction either 'fwd' (initial anchor) or 'back' (final anchor)
% shortfrag coordinates of peudo-torsion defining atoms of the library
%           fragments
%
% codes     vector of allowed fragment indices, empty if none is allowed
%
% G. Jeschke, 25.12.2017

thresh = 0.5; % threshold for coordinate rmsd of pseudotorsion atoms to declare a fragment acceptable

switch direction
    case 'fwd'
        indices = 1:3;
    case 'back'
        indices = 2:4;
end
possible = 0;

codes = zeros(1,length(shortfrag));

for k = 1:length(shortfrag)
    cfrag = shortfrag{k};
    rmsd = superimpose_3points(ptac,cfrag(indices,:));
    if rmsd < thresh
        possible = possible + 1;
        codes(possible) = k;
    end
end

codes = codes(1:possible);
