function [attribute,exceptions] = cx_set_atom(id,address,whichone,attribute)
%
% CX_SET_ATOM Set an attribute of an atom in ChimeraX
%
%   attribute = CX_SET_ATOM(id,address,whichone,attribute)
%   Sets values to an atom atrribute in a ChimeraX entity
%
%   [attribute,exceptions] = CX_SET_ATOM(id,address,whichone,attribute)
%
% INPUT
% id        ChimeraX internal identifier for entity
% address   MMMx address for the atom, see:
%           https://www.cgl.ucsf.edu/chimerax/docs/user/commands/atomspec.html
% whichone  valid attribute name
%           bfactor       double
%           color         https://www.cgl.ucsf.edu/chimerax/docs/user/commands/colornames.html
%                         or (1,3) double RGB triple (0... 255)
%           colorscheme   'byelement', 'byhet'
%           element       string
%           hide          attribute can be missing or empty, otherwise Boolean flag
%           occupancy     double (0 < occupancy <= 1)
%           radius        double
%           selected      Boolean flag
%           show          attribute can be empty, otherwise flag
%           transparency  in %, 0% opaque, 100% fully transparent
% attribute attribute value (string) that should be set
%
% OUTPUT
% attribute  echoes the input attribute if successful, is empty otherwise
% exceptions cell vector of MException objects if something went wrong, 
%            defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% ChimeraX requires plural for levels, but MMMx accepts also singular

% allow for capitalization in attribute name
whichone = lower(whichone);

[spec,exceptions] = cx_from_mmmx_address(address,id);
if ~isempty(exceptions{1})
    return
end

switch whichone
    case 'bfactor'
        attribute = sprintf('%8.3f',attribute);
    case 'color'
        if isnumeric(attribute)
            attribute = sprintf('rgb(%i,%i,%i)',attribute); 
        end
    case 'colorscheme'
        attribute = lower(attribute);
    case 'hide'
        if ~exist('attribute','var') || isempty(attribute)
            attribute = 1;
        end
        attribute = sprintf('%i',attribute);
    case 'occupancy'
        attribute = sprintf('%8.4f',attribute);
    case 'radius'
        attribute = sprintf('%8.4f',attribute);
    case 'selected'
        switch attribute
            case 0
                attribute = 'False';
            case 1
                attribute = 'True';
        end
    case 'show'
        if ~exist('attribute','var') || isempty(attribute)
            attribute = 0;
        end
        attribute = sprintf('%i',~attribute);
        whichone = 'hide';
    case 'transparency'
        attribute = sprintf('%i',attribute); 
   otherwise
        exceptions{1} = MException('cx_set_atom:no_such_attribute', 'Attribute %s not supported for an atom', whichone);
end

if strcmpi(whichone,'color') || strcmpi(whichone,'colorscheme')
    cx_command(sprintf('color %s %s',spec,attribute));
elseif strcmpi(whichone,'transparency')
    cx_command(sprintf('transparency %s %s target a',spec,attribute));
else
    attribute = cx_set_attribute(spec,'atom',whichone,attribute);
end
