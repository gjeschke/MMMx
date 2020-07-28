function [entity,exception,id] = cx_get_pdb(ident,options)
%
% CX_GET_PDB Load PDB structure into ChimeraX and retrieve entity
%
%   response = CX_GET_PDB
%   Returns ChimeraX model number and entity in MMMx:atomic representation
%
%   [entity,exception,id] = CX_GET_PDB(ident,options)
%   Loads structure from PDB server or local file into ChimeraX vi its
%   'open' command and returns the structure as an entity in MMMx:atomic
%   format, for future reference, the ChimeraX identifier is returned
%
% INPUT
% ident     ChimeraX comman, see: 
%           https://www.cgl.ucsf.edu/chimerax/docs/user/index.html
% port      optional port number, defaults to 51051
%           if you override the default, use a value N in the range 
%           49152-65535 and make sure that ChimeraX executes at startup:
%           remotecontrol rest start port N
% options   optional modifications of connection setup, see Matlab
%           documentation of 'weboptions'
%
% OUTPUT
% response  ChimeraX response, a single long string that may contain
%           carriage returns and may be empty, the line after the last
%           carriage return is empty
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

id = [];
entity = [];
exception = [];

if ~exist('ident','var') || isempty(ident)
    fname = 'browse';
end

[path,basname,ext] = fileparts(fname);

if ~isempty(path)
    response = ChimeraX(sprintf('cd %s',path));
    fprintf(1,'%s\n',response);
elseif ~isempty(ext) || (length(fname) ~= 4 && ~strcmp(fname,'browse'))
    response = ChimeraX(sprintf('cd %s',pwd));
    fprintf(1,'%s\n',response);
end


response = ChimeraX(sprintf('open %s%s',basname,ext));
fprintf(1,'%s%\n',response);

poi = strfind(response,'Chain information for');
if isempty(poi)
    return
end
idstr= split(response(poi:end),'#');
idstr = strtok(idstr);
idstr = idstr{2};
poi = strfind(idstr,'.');
if ~isempty(poi)
    idstr = idstr(1:poi-1);
end
id = str2double(idstr);
if isempty(id) || isnan(id)
    return
end

response = ChimeraX(sprintf('getcrd #%i',id));
tic,
atoms = splitlines(response);
for k = 1:length(atoms)
    if ~isempty(atoms{k})
        atom_info = textscan(atoms{k},'Atom /%s:%f@%s%f%f%f');
        atname = atom_info{3};
        entity.(atom_info{1}).(sprintf('R%i',atom_info{2})).(atname{1})=...
            [atom_info{4},atom_info{5},atom_info{6}];
    end
end
toc,

disp('Aber hallo!');
