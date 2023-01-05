function [attribute,exceptions] = cx_set_chain(id,address,whichone,attribute)
%
% CX_SET_CHAIN Set an attribute of a chain in ChimeraX
%
%   attribute = CX_SET_CHAIN(id,address,whichone,attribute)
%   Sets values to a chain atrribute in a ChimeraX entity
%
%   [attribute,exceptions] = CX_SET_CHAIN(id,address,whichone,attribute)
%
% INPUT
% id        ChimeraX internal identifier for entity
% address   MMMx address for the residue
% whichone  valid attribute name
%           color         https://www.cgl.ucsf.edu/chimerax/docs/user/commands/colornames.html
%                         or (1,3) double RGB triple (0... 255)
%           hide          attribute can be missing or empty, otherwise Boolean flag
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
    case 'color'
        if isnumeric(attribute)
            attribute = sprintf('rgb(%i,%i,%i)',attribute); 
        end
    case 'hide'
        if ~attribute
            whichone = 'show';
        end
    case 'selected'
    case 'show'
        if ~attribute
            whichone = 'hide';
        end
    case 'transparency'
        attribute = sprintf('%i',attribute); 
   otherwise
        exceptions{1} = MException('cx_set_chain:no_such_attribute', 'Attribute %s not supported for a chain', whichone);
end

if strcmpi(whichone,'color')
    cx_command(sprintf('color %s %s',spec,attribute));
elseif strcmpi(whichone,'transparency')
    cx_command(sprintf('transparency %s %s target c',spec,attribute));
elseif strcmpi(whichone,'selected')
    if attribute
        cx_command(sprintf('select %s',spec));
    else
        cx_command(sprintf('~select %s',spec));
    end
elseif strcmpi(whichone,'show') || strcmpi(whichone,'hide')
    cx_command(sprintf('%s %s target c',whichone,spec));
end
