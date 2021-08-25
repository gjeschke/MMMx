function complete = check_sidechain_integrity(residue)
%  complete = CHECK_SIDECHAIN_INTEGRITY(residue)
%    Check whether a sidechain contains all expected heavy atoms
%
% INPUT:
%
% residue   MMMx:atomic residue field
%
% OUTPUT:
%
% complete  flag indicating whether the residue features all expected heavy
%           atoms
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

defs = load('monomer_attributes.mat');
aa = defs.monomers.aa;
aa_tags = defs.monomers.aa_tags;

complete = true;

rid = tag2id(upper(residue.name),upper(aa_tags));
if isempty(rid) % non-amino acid residues are declared complete
    return
end
atoms = aa.atoms{rid};
aid = 0;
atag = '?';

while ~isempty(atag)
    aid = aid + 1;
    atag = id2tag(aid,atoms);
    if ~isempty(atag) && atag(1) ~= 'H' && ~strcmp(atag,'OXT') % the remaining tags should exist
        found = false;
        ratoms = fieldnames(residue);
        for ka = 1:length(ratoms) % expand over all atoms
            atom = ratoms{ka};
            if isstrprop(atom(1),'upper') % these are atom fields
                if strcmpi(atom,atag)
                    found = true;
                end
            end
        end
        if ~found
            complete = false;
        end
    end
end