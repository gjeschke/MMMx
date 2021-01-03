function slc = tlc2slc(tlc)
% slc = TLC2SLC(tlc)
% Converts amino acid three-letter code to single-letter code

slc = ''; % initialize empty output

tlcs = upper(':Ala:Arg:Asn:Asp:Cys:Gln:Glu:Gly:His:Ile:Leu:Lys:Met:Phe:Pro:Ser:Thr:Trp:Tyr:Val:Asx:Glx:Xaa:');
slcs = 'ARNDCQEGHILKMFPSTWYVBZX';

id = tag2id(upper(tlc),tlcs);
if isempty(id)
    return
end
slc = slcs(id);

function id=tag2id(tag,tags)
% function id=tag2id(tag,tags)
%
% Returns the identification code corresponding to a string tag by
% comparison with a list of possible tags given in a string with colon (:)
% separation
% if tag is not found in tags or if the position is larger than the number
% of elements in codelist, an empty id is given back
%
% tag       (string) tag to be found, example: 'Be'
% tags      (string) tag list, example ':H:He:Li:Be:C:N:O:F:Ne:'
%
% id        integer telling the position of tag in tags, for example above:
%           id=4
%
% G. Jeschke, 2009

id=[]; % empty default output
delimiter=':'; % colon as default delimiter

etag=[delimiter tag delimiter];
position=strfind(tags,etag);
if position
    id=1+sum(find(tags==delimiter)<position(1));
end

