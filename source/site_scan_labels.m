function [entity,failures] = site_scan_labels(entity,label,fname,options)
%
% SITE_SCAN_LABELS Spin labelling site scan
%
%   [entity,failures] = SITE_SCAN_LABELS(entity,label,fname,options)
%   Finds potential spin labelling sites matching certain conditions and
%   saves them in a list file
%
% INPUT
% entity    entity in an MMMx format
% label     label name, such as 'mtsl' or 'I1A'
% fname     file name for the output list, extension '.lst' is appended if
%           none is provided
% options   computation options
%           .restypes     allowed residue types for labeling, string of
%                         single-letter codes, '*' allows for all amino
%                         acids, defaults to '*'
%           .min_rotamers minimum number of rotamers, defaults to 1 
%           .min_Z        minimum partition function, defaults to 0.1 
%           .chains       restrict to chains, given in a string of
%                         identifiers, '*' denotes all chains, defaults to
%                         all chains
%           .ensemble     if true, loop over all conformers of the entity,
%                         defaults to false (only first conformer)
%
% OUTPUT
% entity       the input entity augmented by newly computed labels
% failures     struct with fields that count failed sites
%              .impossible  sites where no rotamers were returned
%              .too_few     sites with less than min_rotamers rotamers
%              .loo_low_Z   sites with partition function below than min_Z 
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2010-2020: Yevhen Polyhach, Gunnar Jeschke

% initialize output

failures.impossible = 0;
failures.too_few = 0;
failures.too_low_Z = 0;

% initialize options, if it does not exist
if ~exist('options','var') || isempty(options)...
        || ~isfield(options,'restypes') || isempty(options.restypes)
    options.restypes = '*';
end

% set defaults of computation options if required
if ~isfield(options,'min_rotamers') || isempty(options.min_rotamers)
    options.min_rotamers = 1;
end
if ~isfield(options,'min_Z') || isempty(options.min_Z)
    options.min_Z = 0.1;
end
if ~isfield(options,'chains') || isempty(options.chains)
    options.chains = '*';
end
conformers = 1;
if isfield(options,'ensemble') && ~isempty(options.ensemble) && options.ensemble
    conformers = 1:length(entity.populations);
end


% fix output file name
if isempty(fname)
    fname = 'MMMx_sitescan.lst';
else
    [pname,fname_core,ext] = fileparts(fname);
    if isempty(ext)
        fname = fullfile(pname,[fname_core '.lst']);
    else
        fname = fullfile(pname,[fname_core ext]); 
    end
end

fid = fopen(fname,'wt');
fprintf(fid,'%% MMMx spin labeling site scan for chains (%s) and label: %s\n',options.chains,label);
fprintf(fid,'%% Site       rotamers  partition function  position rmsd (%s)    residue\n',char(197));

% scan over all chains 
chains = fieldnames(entity);
for kc = 1:length(chains)
    chain = chains{kc};
    if isstrprop(chain(1),'upper') % chain fields start with a capital
        if strcmp(options.chains,'*') || contains(options.chains,chain(1)) % the chain is selected
            residues = fieldnames(entity.(chain));
            for kr = 1:length(residues)
                residue = residues{kr};
                if strcmp(residue(1),'R') % these are residue fields
                    slc = tlc2slc(entity.(chain).(residue).name);
                    if isempty(slc)
                        slc = na2slc(entity.(chain).(residue).name);
                    end
                    if ~isempty(slc) % otherwise, it is not an amino acid or nucleic acid
                        if strcmp(options.restypes,'*') || contains(options.restypes,slc) % the amino acid is selected
                            for c = conformers
                                site = sprintf('{%i}(%s)%s',c,chain,residue(2:end));
                                [argsout,entity,exceptions] = get_label(entity,label,{'positions','populations','part_fun'},site);
                                if ~isempty(exceptions{1})
                                    failures.impossible = failures.impossible + 1;
                                    continue
                                end
                                positions = argsout{1}{1};
                                populations = argsout{2}{1};
                                nrot = length(populations);
                                if nrot < options.min_rotamers
                                    failures.too_few = failures.too_few + 1;
                                    continue;
                                end
                                Z = argsout{3}{1};
                                if Z < options.min_Z
                                    failures.too_low_Z = failures.too_low_Z + 1;
                                    continue;
                                end
                                populations = populations/sum(populations);
                                mean_coor = populations'*positions;
                                coor_dev = positions - repmat(mean_coor,nrot,1);
                                rmsd = sqrt(sum(sum(coor_dev.^2))/nrot);
                                nrot_str = pad(sprintf('%i',nrot),21-length(site),'left');
                                res_str = pad(sprintf('{%s}',entity.(chain).(residue).name),5,'left');
                                fprintf(fid,'%s%s%20.3f%19.2f      %s\n',site,nrot_str,Z,rmsd,res_str);
                            end
                        end
                    end
                end
            end
        end
    end
end

fclose(fid);