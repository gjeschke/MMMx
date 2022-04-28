function [entity,failures] = hetero_pair_list(entity,sites1,label1,sites2,label2,fname,options)
%
% HETERO_PAIR_LIST Label pair list for different labels
%
%   [entity,failures] = HETERO_PAIR_LIST(entity,sites1,label1,sites2,label2,fname,options)
%   Finds potential spin labelling sites matching certain conditions and
%   saves them in a list file
%
% INPUT
% entity    entity in an MMMx format
% sites1    structure with tested labelling sites, only field .address is
%           used
% label1    label name, such as 'mtsl' or 'I1A' for sites1
% sites2    structure with tested labelling sites, only field .address is
%           used
% label2    label name, such as 'mtsl' or 'I1A' for sites2
% fname     file name for the output list, extension '.lst' is appended if
%           none is provided
% options   computation options
%           .r_min        minimum mean distance, defaults to 15 Angstroem 
%           .r_max        maximum mean distance, defaults to 80 Angstroem 
%           .coupled      if true, only pairs within the same conformer are
%                         considered, defaults to true
%
% OUTPUT
% entity       entity in an MMMx format, may contain new label rotamers
% failures     struct with fields that count failed sites
%              .empty       pairs has an empty distance distribution within
%                           the sensitive range
%              .too_short   pairs with too short mean distance
%              .too_long    pairs with too long mean distance
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2010-2020: Yevhen Polyhach, Gunnar Jeschke

% initialize output

failures.empty = 0;
failures.too_short = 0;
failures.too_long = 0;

% initialize options, if it does not exist
if ~exist('options','var') || isempty(options)...
        || ~isfield(options,'r_min') || isempty(options.r_min)
    options.r_min = 15;
end

% set defaults of computation options if required
if ~isfield(options,'r_max') || isempty(options.r_max)
    options.r_max = 80;
end
if ~isfield(options,'coupled') || isempty(options.coupled)
    options.coupled = true;
end


% fix output file name
if isempty(fname)
    fname = 'MMMx_pairlist.lst';
else
    [pname,fname_core,ext] = fileparts(fname);
    if isempty(ext)
        fname = fullfile(pname,[fname_core '.lst']);
    else
        fname = fullfile(pname,[fname_core ext]); 
    end
end

fid = fopen(fname,'wt');
fprintf(fid,'%% MMMx site pair list for labels: %s|%s\n',label1,label2);
fprintf(fid,'%% Site1            Site2     rmean (%s)  sigr (%s)\n',char(197),char(197));

S1 = length(sites1);
S2 = length(sites2);

for s1 = 1:S1-1
    c1 = which_conformer(sites1(s1).address);
    for s2 = 1:S2
        c2 = which_conformer(sites2(s2).address);
        if c1 == c2 || ~options.coupled
            [r_axis,distribution,entity] = distance_distribution(entity,...
                sites1(s1).address,label1,sites2(s2).address,label2,options);
            if isempty(distribution)
                failures.empty = failures.empty + 1;
                continue;
            end
            distribution = distribution/sum(distribution);
            rmean = sum(r_axis.*distribution);
            if rmean < options.r_min
                failures.too_short = failures.too_short + 1;
                continue;
            end
            if rmean > options.r_max
                failures.too_long = failures.too_long + 1;
                continue;
            end
            Delta_r = r_axis - rmean;
            sigr = sqrt(sum(distribution.*Delta_r.^2));
            site1 = sites1(s1).address;
            site2_str = pad(sites2(s2).address,24-length(site1),'left');
            fprintf(fid,'%s%s%14.1f%10.1f\n',site1,site2_str,rmean,sigr);
        end
    end
end

fclose(fid);

function c = which_conformer(address)
% retrieves conformer number from site address

c = 1;
pa = strfind(address,'{');
pe = strfind(address,'}');
if ~isempty(pa) && ~isempty(pe)
    c = str2double(address(pa+1:pe-1));
end