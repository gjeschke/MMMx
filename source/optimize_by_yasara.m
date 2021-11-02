function exceptions = optimize_by_yasara(logfile,fname,options)
% function optimize_by_yasara(logfile,fname,options)
%
% Optimizes the structure of proteins, nucleic cids, and their complexes by
% calling Yasara, which must be on the Matlab path
%
% Input:
%
% logfile   file identifier for log file
% fname     file name for input file, output file name is input file name 
%           with '_yasara' appended, unless specified in options
% options   structure with options specified by fields
%           .console    flag, if true, Yasara is run in console mode
%           .fname      if specified, output file name
%
%
% Output:
%
% exceptions    cell(1,1) of a potential Matlab exception, if input or
%               output file could not be opened
%
% G. Jeschke, 7.5.2021

exceptions{1} = [];

% Yasara console mode defaults to false
if ~exist('options','var') || isempty(options) || ~isfield(options,'console')
    options.console = false;
end

% if no maximum time is specified, the Yasara job is limited to 1 hour
if ~isfield(options,'max_time') || isempty(options.max_time)
    options.max_time = 1;
end

% find Yasara executable
yasara = which_third_party_module('yasara');
if isempty(yasara)
    fprintf(logfile,'ERROR: Yasara is not on path. Skipping optimization');
    return
end

% find directory where Yasara is located
yasara_path = fileparts(yasara); 

% analyze file name and add current directory, if required
[fpath,bname,ext] = fileparts(fname);
if isempty(ext)
    ext = '.pdb';
end
if isempty(fpath)
    fpath = pwd;
end

% specify full path name of input file
ffname = fullfile(fpath,strcat(bname,ext));

% determine output file name
if isfield(options,'fname') && ~isempty(options.fname)
    yname = options.fname;
    if ~contains(yname,'.pdb')
        yname = strcat(yname,'.pdb');
    end
else
    yname = strcat(bname,'_yasara.pdb');
end
% specify full path name of output file
foname = fullfile(fpath,yname);

mypath = pwd; % store current path
cd(yasara_path); % change to Yasara directory

% read the Yasara optimization template control file and adapt it to
% current file name
ctrl_file = sprintf('%s.mcr',bname);
iname = 'minimization_server.mcr';
if ~isdeployed
    iname = which('minimization_server.mcr');
end
ifid = fopen(iname,'r');
ofid = fopen(ctrl_file,'wt');
fprintf(ofid,'LoadPDB %s\n',ffname);
while 1
    tline = fgetl(ifid);
    if ~ischar(tline), break, end
    fprintf(ofid,'%s\n',tline);
end
% fprintf(ofid,'SelectRes %i-%i\n',res(1),res(2));
fprintf(ofid,'SavePDB 2,%s_opt.pdb,Format=PDB,Transform=Yes\n',bname);
fprintf(ofid,'Exit\n');
fclose(ifid);
fclose(ofid);
cmd = sprintf('yasara %s &',ctrl_file);
if isunix
    cmd = strcat('./',cmd);
end
if options.console
    cmd = strcat(cmd,' -con');
end
cmd = strcat(cmd,' &');

y_output = sprintf('%s_opt.pdb',bname);
if exist(y_output,'file')
    delete(y_output);
end
runclock = tic;
t = timer;
t.ExecutionMode = 'fixedSpacing';
t.Period = 5;
t.StartDelay = 5;
t.UserData.y_output = y_output;
t.UserData.runclock = runclock;
t.UserData.maxtime = 3600*options.max_time;
tasks = ceil(t.UserData.maxtime/t.Period);
t.TasksToExecute = tasks;
t.UserData.status = 'running';
t.TimerFcn = @check_for_completion;
start(t);
[s, w] = system(cmd);
if s~=0
    fprintf(logfile,'\nERROR: Yasara did not complete successfully:\n%s\n\n',w);
    delete(ctrl_file);
    exceptions = {MException('optimize_by_yasara:yasara_failed',...
                                    'Yasara run failed (see logfile)')};
    cd(mypath);
    return
end
wait(t);
fprintf(logfile,'Yasara job %s\n',t.UserData.status);
delete(t);

% the following is required, because Yasara does not produce PDB standard
% format
if exist(y_output,'file')
    exceptions = correct_yasara(y_output,foname);
    delete(y_output);
end
delete(ctrl_file);


% change back to original directory
cd(mypath);

function check_for_completion(mTimer,~)

mytime = toc(mTimer.UserData.runclock);
if exist(mTimer.UserData.y_output,'file')
    mTimer.UserData.status = sprintf('completed in %4.1f s.',mytime);
    stop(mTimer);
else    
    if mytime > mTimer.UserData.maxtime
        mTimer.UserData.status = 'timeout';
        stop(mTimer);
    end
end

