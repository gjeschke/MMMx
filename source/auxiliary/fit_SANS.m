function [fit,chi2,outname,status,result] = fit_SANS(datafile,pdbfile,illres,options)
%
% FIT_SANS Fits a small-angle neutron scattering curve by calling CRYSON
%          from the ATSAS package
%
%   [fit,chi2,outname,status,result] = FIT_SANS(datafile,pdbfile,illres,options)
%          requires CRYSON to be on the Matlab path, returns
%          empty output if not
%
% INPUT
% datafile  SANS data file
% pdbfile   PDB file of the structural model
% illres    name of the resolution file, can be missing, but this is not
%           recommended
% options   additional options
%           .lm           maximum order of harmonics (1...50), default: 15
%           .fb           order of Fibonacci grid (10...18), default: 17
%           .D2O          heavy water content in the solution 0<= .D2O<= 1,
%                         default: 0 
%           .per          Boolean flag, if true, chain is assumed to be
%                         perdeuterated, default: false
%           .err          Bolean flag, if true, error is returned as last
%                         column with sufficiently new versions
%           .un           angular units of experimental data, defaults to
%                         automatic estimate
%                         1     4 pi sin(theta)/Angstroem
%                         2     4 pi sin(theta)/nm
%                         3     2 sin(theta)/Angstroem
%                         4     2 sin(theta)/nm
%           .cst          Boolean flag, if true constant subtraction for
%                         mismatched buffer is done, defaults to true
%           .eh           Boolean flag, if true account for explicit
%                         hydrogens
%           .delete       Boolean flag, if tru, delete output file after
%                         reading fit, default: true
%
% OUTPUT
% fit        double (ndat,3) or (ndat,4) with fit result
%            column 1: angular axis
%            column 2: experimental data
%            column 3: fit curve
%            column 4: error of experimental data
% chi2       chi^2 value
% outname    defaults to one cell holding an empty array
% status     status answer of the CRYSON call
% result     result answer of the CRYSON call
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

if ~exist('options','var')
    options.err = false;
end

poi = strfind(pdbfile,'.pdb');
if isempty(poi)
    outname = strcat(pdbfile,'00.fit');
    pdbfile = strcat(pdbfile,'.pdb');
    to_be_deleted = strcat(pdbfile,'00*.*');
else
    outname = [pdbfile(1:poi-1) '00.fit'];
    to_be_deleted = strcat(pdbfile(1:poi-1),'00*.*');
end

s=which('cryson.exe');

if isempty(s)
    return
end

cmd=[s ' ' pdbfile ' ' datafile];
if exist('illres','var') && ~isempty(illres)
    cmd = [cmd ' -res ' illres];
end

if isfield(options,'lm')
    cmd = sprintf('%s -lm %i',cmd,options.lm);
end

if isfield(options,'fb')
    cmd = sprintf('%s -fb %i',cmd,options.fb);
end

if isfield(options,'D2O')
    cmd = sprintf('%s -D2O %5.2f',cmd,options.D2O);
end

if isfield(options,'per') && options.per
    cmd = sprintf('%s -per',cmd);
end

if isfield(options,'err') && options.err
    cmd = sprintf('%s -err',cmd);
end

if isfield(options,'un') && ~isempty(options.un)
    cmd = sprintf('%s -un %i',cmd,options.un);
end

if isfield(options,'cst') && ~options.cst
    cmd = sprintf('%s -cst',cmd);
end

if isfield(options,'eh') && ~options.eh
    cmd = sprintf('%s -eh',cmd);
end

[status,result]=dos(cmd);

fit = zeros(10000,4);

chi2 = [];

fid = fopen(outname);
if fid==-1
    fit = [];
    return;
end
nl=0;
while 1
    tline = fgetl(fid);
    if ~ischar(tline), break, end
    if nl > 0 % skip first line
        dataset = str2num(tline); %#ok<ST2NM>
        ncol = length(dataset);
        fit(nl,1:ncol) = dataset;
    else
        poi = strfind(tline,'Chi^2:');
        if ~isempty(poi)
            rem = textscan(tline(poi+6:end),'%s');
            args = rem{1};
            chi2 = str2double(char(args(1)));
        else
            poi = strfind(tline,'Chi:');
            rem = textscan(tline(poi+4:end),'%s');
            args = rem{1};
            chi = str2double(char(args(1)));
            chi2 = chi^2;
        end
    end
    nl=nl+1;
end
fit = fit(1:nl-1,:);
if ncol == 4 % shift error column to end in output of newer CRYSON or remove it
    if options.err
        fit = fit(:,[1 2 4 3]);
    else
        fit = fit(:,[1 2 4]);
    end
end
fclose(fid);

if isfield(options.delete) && options.delete
    delete(to_be_deleted);
end