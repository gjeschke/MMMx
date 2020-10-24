function [attribute,exceptions] = cx_get_atom(id,address,whichone)
%
% CX_GET_ATOM Set an attribute of an atom in ChimeraX
%
%   attribute = CX_GET_ATOM(id,address,whichone,attribute)
%   Gets value of an atom atrribute in a ChimeraX entity
%
%   [attribute,exceptions] = CX_GET_ATOM(id,address,whichone,attribute)
%
% INPUT
% id        ChimeraX internal identifier for entity
% address   MMMx address for the atom
% whichone  valid attribute name
%           alt_locs      cell of location tags
%           bfactor       double
%           color         (1,3) int RGB 0... 255
%           coord         (1,3) double
%           element       string, uppercase
%           hide          flag
%           occupancy     double
%           radius        double
%           selected      flag
%           transparency  double
%
% OUTPUT
% attribute  attribute value, empty if the attribute is not supported or
%            cannot be retrieved
% exceptions cell vector of MException objects if something went wrong, 
%            defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% ChimeraX requires plural for levels, but MMMx accepts also singular

% initialize empty exception cell vector

% allow for capitalization in attribute name
whichone = lower(whichone);

[spec,exceptions] = cx_from_mmmx_address(address,id);
if ~isempty(exceptions{1})
    return
end

% treat syntax difference between MMMx and ChimeraX
if strcmp(whichone,'coor')
    whichone = 'coord';
end

% special treatment fro transparency
if strcmp(whichone,'transparency')
    attribute = cx_get_attribute(spec,'atom','color');
else
    attribute = cx_get_attribute(spec,'atom',whichone);
end



switch whichone
    case 'alt_locs'
        attribute = split(attribute,',');
    case 'bfactor'
        attribute = str2double(attribute);
    case 'color'
        attribute = str2num(attribute); %#ok<ST2NM>
        attribute = attribute(1:3);
    case 'coord'
        attribute = str2num(attribute); %#ok<ST2NM>
    case 'element'
        attribute = upper(attribute);
    case 'hide'
        attribute = str2double(attribute);
    case 'occupancy'
        attribute = str2double(attribute);
    case 'radius'
        attribute = str2double(attribute);
    case 'selected'
        switch attribute
            case 'False'
                attribute = 0;
            case 'True'
                attribute = 1;
        end
    case 'transparency'
        attribute = str2num(attribute); %#ok<ST2NM>
        attribute = 100 - 100*attribute(4)/255;
    otherwise
        exceptions{1} = MException('cx_get_atom:no_such_attribute', 'Attribute %s not supported for an atom', whichone);
end