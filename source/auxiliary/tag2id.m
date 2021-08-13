function id = tag2id(tag,tags,codelist,delimiter)
% id = TAG2ID(tag,tags,codelist,delimiter)
%
% Returns the identification code corresponding to a string tag by
% comparison with a list of possible tags given in a string with colon (:)
% separation, another separator can be selected by optional parameter
% delimiter
% if an array codelist is given the identifier is of the same type as the
% elements in codelist, the sequence in the codelist must correspond to the
% one in the tags string
% if codelist is missing, the id is an integer corresponding to the
% position of the tag in string tags
% if tag is not found in tags or if the position is larger than the number
% of elements in codelist, an empty id is given back
%
% tag       (string) tag to be found, example: 'Be'
% tags      (string) tag list, example ':H:He:Li:Be:C:N:O:F:Ne:'
% codelist  [optional] (array of ids), example [1,4,7,9,12,14,16,19,20],
%           can be empty to allow different delimiter without codelist
% delimiter [optional] (char), alternative delimiter, defaults to ':'
%
% id        id selected from array codelist or integer telling the position
%           of tag in tags, for example above:
%           if codelist is missing, id=4
%           if codelist is present, id=9
%
% G. Jeschke, 2009

id=[]; % empty default output
if nargin<4,
    delimiter=':'; % colon as default delimiter
end;

etag=[delimiter tag delimiter];
position=strfind(tags,etag);
if position,
    id=1+sum(find(tags==delimiter)<position(1));
    if nargin>=3 && ~isempty(codelist),
        if id<=length(codelist),
            id=codelist(id);
        else
            id=[];
        end;
    end;
end;

