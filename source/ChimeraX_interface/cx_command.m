function response = cx_command(command,port,options)
%
% CX_COMMAND Send command to ChimeraX and receive response
%
%   response = CX_COMMAND
%   Returns an (entity) structure in MMMx:atomic representation
%
%   [entity,exception] = CX_COMMAND(command,port,options)
%   Sends a command to the REST interface of ChimeraX and returns the
%   response of ChimeraX in a string variable
%
% INPUT
% command   ChimeraX command, see: 
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


if exist('port','var') && ~isempty(port)
    IP_port = sprintf('http://127.0.0.1:%i',port);
else % default port number
    IP_port = 'http://127.0.0.1:51051';
end

url = sprintf('%s/run?command=%s',IP_port,urlencode(command));

if exist('options','var')
    response = webread(url,options);
else
    response = webread(url);
end