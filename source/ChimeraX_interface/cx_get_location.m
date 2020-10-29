function [attribute,exceptions] = cx_get_location(id,address,whichone)
%
% CX_GET_LOCATION Set an attribute of an object to ChimeraX
%
%   attribute = CX_GET_LOCATION(id,address,whichone,attribute)
%   Sets values to a location atrribute in a ChimeraX entity
%
%   [attribute,exceptions] = CX_GET_LOCATION(id,address,whichone,attribute)
%
% INPUT
% id        ChimeraX internal identifier for entity
% address   MMMx address for the location
% whichone  valid attribute name, only 'coor' is supported for a location
%
% OUTPUT
% attribute  (1,3) double coordinates (?)
% exceptions cell vector of MException objects if something went wrong, 
%            defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% ChimeraX requires plural for levels, but MMMx accepts also singular

attribute = '';

% initialize empty exception cell vector
exceptions = {[]};

% separate location tag
poi = strfind(address,':');
if isempty(poi) || length(address) == poi
    exceptions{1} = MException('cx_get_location:no_location_tag', 'Location tag missing in address: %s', address);
    return
end
location = address(poi+1:end);
address = address(1:poi-1);

[spec,exceptions] = cx_from_mmmx_address(address,id);
if ~isempty(exceptions{1})
    return
end

response = cx_set_attribute(spec,'atom','alt_loc',location);
if isempty(response)
    exceptions{1} = MException('cx_get_location:no_such_location', 'Location %s does not exist for atom %s', location, address);
end

switch whichone
    case 'coor'
        attribute = cx_get_attribute(spec,'atom','coord');
        attribute = str2num(attribute); %#ok<ST2NM>
    otherwise
        exceptions{1} = MException('cx_get_location:no_such_attribute', 'Attribute %s not supported for a location', whichone);
end