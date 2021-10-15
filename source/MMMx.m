function MMMx(control_file)
%
% MMMX Run an MMMx modelling job
%
%   MMMx(control_file)
%   Runs the MMMx job specified by a control file or opens a window for
%   selecting such a control file and runs the job selected by the user
%
% INPUT
% control_file control file in MMMx format, if missing, and 'open file'
%              window is opened
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

my_path = pwd;

if ~exist('control_file','var') || isempty(control_file)
    [control_file,pathname] = uigetfile('*.mcx');
    cd(pathname);
elseif strcmpi(control_file,'help') || strcmpi(control_file,'docs')
    web('index.html');
    return
end

exceptions = MMModel(control_file);

for k = 1:length(exceptions)
    if ~isempty(exceptions{k})
        c_exception = exceptions{k};
        fprintf(2,'%s\n',c_exception.message);
    end
end

cd(my_path);