function [spec,exceptions] = cx_from_mmmx_address(address,id)
%
% CX_FROM_MMMX_ADDRESS Translate ChimeraX target specification
%
%   [spec,exceptions] = CX_FROM_MMMX_ADDRESS(address)
%   Returns ChimeraX address
%
%   [addresses,exceptions] = CX_FROM_MMMX_address(address,id)
%   Projects an MMMx address onto a ChimeraX target specification,
%   and assigns a structure ID
%
% INPUT
% address    MMMx address referring to one entity
% id         optional structure identifier (integer) in ChimeraX, if
%            missing, all structures loaded into ChimeraX are included
%
% OUTPUT
% spec       ChimeraX target specification (string), see:
%            https://www.cgl.ucsf.edu/chimerax/docs/user/commands/atomspec.html
% exceptions cell vector of MException objects if something went wrong, 
%            defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initalize empty output
spec = '';
exceptions = {[]};
warnings = 0; % initialize counter for warnings

if ~exist('id','var')
    id = [];
end

if ~isempty(id)
    spec = sprintf('#%i',id);
end

% extract and encode conformer specification, if any
si = strfind(address,'{');
ei = strfind(address,'}');
if ~isempty(si) && ~isempty(ei)
    conf_str = address(si+1:ei-1);
    if isempty(spec)
        spec = '#*';
    end
    spec = sprintf('%s.%s',spec,conf_str);
    % remove conformer part
    address = [address(1:si-1) address(ei+1:end)];
end

% extract and encode chain specification, if any
si = strfind(address,'(');
ei = strfind(address,')');
if ~isempty(si) && ~isempty(ei)
    chain_str = address(si+1:ei-1);
    spec = sprintf('%s/%s',spec,chain_str);
    address = address(ei+1:end);
end

% extract and encode residue specification, if any
ei = strfind(address,'.');
if ~isempty(ei)
    residue_part = address(1:ei-1);
    address = address(ei+1:end);
else
    ei = strfind(address,':');
    if ~isempty(ei)
        residue_part = address(1:ei-1);
        address = address(ei+1:end);      
    else
        residue_part = address;
        address = '';
    end
end

% check for (disallowed) rotamer specificatoon
ei = strfind(residue_part,'|');
if ~isempty(ei)
    warnings = warnings + 1;
    exceptions{warnings} = MException('cx_from_mmmx_address:invalid_rotamer_spec', 'Rotamer specification %s is disregarded', residue_part(ei+1:end));
    residue_part = residue_part(1:ei-1);
end

% add residue specification
if ~isempty(residue_part)
    spec = sprintf('%s:%s',spec,residue_part);
end

atom_part = address;

% check for (disallowed) location specificatoon
ei = strfind(atom_part,':');
if ~isempty(ei)
    warnings = warnings + 1;
    exceptions{warnings} = MException('cx_from_mmmx_address:invalid_location_spec', 'Location specification %s is disregarded', atom_part(ei+1:end));
    atom_part = atom_part(1:ei-1);
end

% add atom specification
if ~isempty(atom_part)
    spec = sprintf('%s@%s',spec,atom_part);
end
