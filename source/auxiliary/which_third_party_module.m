function executable = which_third_party_module(module)
% WHICH_THIRD_PARTY_MODULE(module)
%
% finds the pathname of the executable for a third-party module, using the
% registration list third-part.lst
%
% INPUT:
%
% module    name of the requested third-party module, mandatory
%
% OUTPUT:
%
% executable    pathname to the registered executable, empty, if a module
%               of that name is not registered

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% initialize empty output
executable = '';

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

found = false;
for k = 1:length(modules)
    if strcmpi(modules(k).name,module)
        for k2 = 1:length(modules(k).executable)
            [status, result] = system(sprintf('where %s',modules(k).executable{k2}));
            if status == 0
                executable = strtrim(result);
                break
            end
        end
        if ~isempty(executable)
            found = true;
        elseif ~isdeployed
            exe = 0;
            while ~found && exe < length(modules(k).executable)
                exe = exe + 1;
                executable = which(modules(k).executable{exe});
                if ~isempty(executable)
                    found = true;
                end
            end
        end
        if found
            break
        end
    end
end

