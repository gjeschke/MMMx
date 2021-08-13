function tag = id2tag(id,tags,codelist,delimiter)
% tag = ID2TAG(id,tags,codelist,delimiter)
%
% Returns the string tag corresponding to an identification code from
% a list of possible tags given in a string with colon (:)
% separation, another separator can be selected by optional parameter
% delimiter
% if an array codelist is given the identifier is of the same type as the
% elements in codelist, the sequence in the codelist must correspond to the
% one in the tags string
% if codelist is missing, the id is an integer corresponding to the
% position of the tag in string tags
% if the id is larger than the number of tags in the tags string an empty
% tag is given back
% if the id appears twice or more in the code list, the first tag is given
% back
%
% id        identifier, for example 7
% tags      (string) tag list, example ':H:He:Li:Be:C:N:O:F:Ne:'
% codelist  [optional] (array of ids), example [1,4,7,9,12,14,16,19,20]
% delimiter [optional] (char), alternative delimiter, defaults to ':'
%
% tag       (string) tag
%           if codelist is missing, tag for the example will be 'O'
%           if codelist is present, tag for the example will be 'Li'
%
% G. Jeschke, 2009

tag=''; % empty default output
if nargin<4,
    delimiter=':'; % colon as default delimiter
end;

if nargin>2,
    poi=find(id==codelist);
    if ~isempty(poi),
        id=poi(end);
    else
        id=length(tags)+1;
    end;
end;

lookup=find(tags==delimiter);
if id<length(lookup),
    tag=tags(lookup(id)+1:lookup(id+1)-1);
end;
