function attribute = cx_get_attribute(spec,level,whichone)
%
% CX_GET_ATTRIBUTE Get an attribute of an object from ChimeraX
%
%   attribute = CX_GET_ATTRIBUTE(spec,level,whichone)
%   Receives an attribute of an atom, residue, chain, or model
%   from ChimeraX as a string
%
% INPUT
% spec      ChimeraX target specification, see:
%           https://www.cgl.ucsf.edu/chimerax/docs/user/commands/atomspec.html
% level     can be atom, residue, chain, model or others listed at
%           https://www.cgl.ucsf.edu/chimerax/docs/user/commands/info.html
% whichone  any valid attribute name at this level
%
% OUTPUT
% attribute  attribute value(s) as a string, client must interpret, is
%            empty if attribute does not exist
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% ChimeraX requires plural for levels, but MMMx accepts also singular
switch level
    case 'atom'
        level = 'atoms';
    case 'residue'
        level = 'residues';
    case 'chain'
        level = 'chains';
    case 'model'
        level = 'models';
end

attribute = cx_command(sprintf('info %s %s attribute %s',level,spec,whichone));

response = strfind(attribute,whichone);

% extract the actual response
if ~isempty(response)
    attribute = strip(attribute(response+length(whichone):end));
else
    attribute = '';
end