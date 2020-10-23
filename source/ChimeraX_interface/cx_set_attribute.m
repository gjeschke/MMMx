function attribute = cx_set_attribute(spec,level,whichone,attribute)
%
% CX_SET_ATTRIBUTE Set an attribute of an object to ChimeraX
%
%   attribute = CX_SET_ATTRIBUTE(spec,level,attribute)
%   Receives an attribute of an atom, residue, chain, or model
%   from ChimeraX as a string
%
% INPUT
% spec      ChimeraX target specification, see:
%           https://www.cgl.ucsf.edu/chimerax/docs/user/commands/atomspec.html
% level     can be atom, residue, chain, model or others listed at
%           https://www.cgl.ucsf.edu/chimerax/docs/user/commands/info.html
% whichone  any valid attribute name at this level
% attribute attribute value (string) that should be set
%
% OUTPUT
% attribute  echoes the input attribute if successful, is empty otherwise
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

response = cx_command(sprintf('setattr %s %s %s %s',spec,level,whichone,attribute));

% check whether the attribute setting succeeded
if ~contains(response,whichone) || contains(response,'Cannot set attribute')
    attribute = '';
end