function find_third_party_modules(pathname)
% FIND_THIRD_PARTY_MODULES(path)
%
% compiles a list of available third-party executables and stores it as
% third-party.lst in a directory specied by input 'path' or in the current
% directory
%
% INPUT:
%
% path  path to the directory where the list is to be stored, defaults 

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

if ~exist('pathname','var') || isempty(pathname)
    pathname = pwd;
end

fname = fullfile(pathname,'third-party.lst');

modules(1).name = 'crysol';
modules(1).executable{1} = 'crysol.exe';
modules(2).name = 'crysol3';
modules(2).executable{1} = 'crysol_30.exe';
modules(3).name = 'cryson';
modules(3).executable{1} = 'cryson.exe';
modules(4).name = 'dssp';
modules(4).executable{1} = 'dssp-2.0.4-win32.exe';
modules(4).executable{2} = 'dssp.exe';
modules(4).executable{3} = 'dsspcmbi.exe';
modules(5).name = 'hydropro10';
modules(5).executable{1} = 'hydropro10-msd.exe';
modules(6).name = 'msms';
modules(6).executable{1} = 'msms.exe';
modules(7).name = 'muscle';
modules(7).executable{1} = 'muscle.exe';
modules(8).name = 'scwrl4';
modules(8).executable{1} = 'scwrl4.exe';
modules(9).name = 'yasara';
modules(9).executable{1} = 'yasara.exe';

fid = fopen(fname,'wt');
for m = 1:length(modules)
    found = false;
    exe = 0;
    while ~found && exe < length(modules(m).executable)
        exe = exe + 1;
        executable = which(modules(m).executable{exe});
        if ~isempty(executable)
            found = true;
        end
    end
    if found
        fprintf(fid,'%s> %s\n',modules(m).name,executable);
        fprintf(1,'Third-party module %s was found.\n',modules(m).name);
    else
        fprintf(2,'Third-party module %s is not present.\n',modules(m).name);
    end
end
fclose(fid);
