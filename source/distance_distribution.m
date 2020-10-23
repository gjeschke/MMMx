function [r_axis,distribution,entity,exceptions] = distance_distribution(entity,site1,label1,site2,label2,options)
%
% DISTANCE_DISTRIBUTION Distance distribution between labels or atoms
%
%   [r_axis,distribution,entity,exceptions] =
%   DISTANCE_DISTRIBUTION(entity,site1,label1,site2,label2)
%   Distance distribution between residues site1 and site2 with label/atom
%   specifications label1 and label2
%
%   [entity,exceptions] =
%   DISTANCE_DISTRIBUTION(entity,site1,label1,site2,label2,options)
%   Allows to specify options for the computation
%
% INPUT
% entity    entity in an MMMx format
% site1     residue address string for site 1 with optional conformer 
%           specifier, use '{*}' for all conformers, defaults to '{1}'
% label1    label name, such as 'mtsl' or 'I1A' or atom specifier for this
%           residue, such as 'atom.CA', for label at site 1
%           specifier, use '{*}' for all conformers, defaults to '{1}'
% site2     residue address string for site 2 with optional conformer 
%           specifier, use '{*}' for all conformers, defaults to '{1}'
% label2    label name, such as 'mtsl' or 'I1A' or atom specifier for this
%           residue, such as 'atom.CA', for label at site 2
% options   computation options
%           .rmin         minimum distance (Angstroem), default 10, double
%           .rmax         maximum distance (Angstroem), default 150, double
%           .resolution   resolution (Angstroem), defaul 0.5, double
%           .units        string, 'probability' makes sum(distribution) = 1
%                         'density' returns probability density in units of
%                         (1/Angstroem), default: 'density'
%           .coupled      boolean, if true (default), only distributions
%                         within the same conformer are added, otherwise
%                         sites in different conformers are also combined
%           .smoothing    standard deviation for Gaussian smoothing,
%                         defaults to twice options.resolution
%
% OUTPUT
% r_axis       distance axis in units of Angstroem
% distribution distance distribution
% entity       the input entity augmented by newly computed labels
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2010-2020: Yevhen Polyhach, Gunnar Jeschke

% initialize output
distribution = [];

exceptions{100} = [];
warnings = 0; % counter for warnings

% initialize options, if it does not exist
if ~exist('options','var') || isempty(options)...
        || ~isfield(options,'rmin') || isempty(options.rmin)
    options.rmin = 10;
end

% set defaults of computation options if required
if ~isfield(options,'rmax') || isempty(options.rmax)
    options.rmax = 150;
end
if ~isfield(options,'resolution') || isempty(options.resolution)
    options.resolution = 0.5;
end
if ~isfield(options,'units') || isempty(options.units)
    options.units = 'density';
end
if ~isfield(options,'coupled') || isempty(options.coupled)
    options.coupled = true;
end
if ~isfield(options,'smoothing') || isempty(options.smoothing)
    options.smoothing = 2*options.resolution;
end

r_axis = options.rmin:options.resolution:options.rmax;

% extract conformer selection of site1 if any, set default, if none
si = strfind(site1,'{');
ei = strfind(site1,'}');
if ~isempty(si) && ~isempty(ei)
    confstr = site1(si+1:ei-1);
    residue1 = [site1(1:si-1) site1(ei+1:end)];
else
    confstr = '1';
end
if strcmp(strtrim(confstr),'*')
    conformers1 =  1:length(entity.populations);
    confstr = '';
end
% make conformer vector
if ~isempty(confstr)
    conformers1 = zeros(1,66000);
    delimited = strfind(confstr,',');
    delimited = [0 delimited length(confstr)+1];
    cpoi = 0;
    for k = 1:length(delimited)-1
        range = confstr(delimited(k)+1:delimited(k+1)-1);
        hyphen = strfind(range,'-');
        if length(hyphen) > 1
            exceptions = {MException('distance_distribution:wrong_conformer_range', 'Conformer range for site 1 contains more than one hyphen')};
            return
        end
        if isempty(hyphen) % single conformer
            cpoi = cpoi + 1;
            % process start/end syntax
            if strcmp(strtrim(range),'start')
                conformers1(cpoi) = 1;
            elseif strcmp(strtrim(range),'end')
                conformers1(cpoi) = length(entity.populations);
            else
                conformers1(cpoi) = str2double(range);
            end
        else % conformer range
            st_str = range(1:hyphen-1);
            % process start/end syntax
            if strcmp(strtrim(st_str),'start')
                st_val = 1;
            elseif strcmp(strtrim(st_str),'end')
                st_val = length(entity.populations);
            else
                st_val = str2double(st_str);
            end
            en_str = range(hyphen+1:end);
            % process start/end syntax
            if strcmp(strtrim(en_str),'start')
                en_val = 1;
            elseif strcmp(strtrim(en_str),'end')
                en_val = length(entity.populations);
            else
                en_val = str2double(en_str);
            end
            crange = st_val:en_val;
            conformers1(cpoi+1:cpoi+length(crange)) = crange;
            cpoi = cpoi + length(crange);
        end
    end
    conformers1 = conformers1(1:cpoi);
end

% extract conformer selection of site1 if any, set default, if none
si = strfind(site2,'{');
ei = strfind(site2,'}');
if ~isempty(si) && ~isempty(ei)
    confstr = site2(si+1:ei-1);
    residue2 = [site2(1:si-1) site2(ei+1:end)];
else
    confstr = '1';
end
if strcmp(strtrim(confstr),'*')
    conformers2 =  1:length(entity.populations);
    confstr = '';
end
% make conformer vector
if ~isempty(confstr)
    conformers2 = zeros(1,66000);
    delimited = strfind(confstr,',');
    delimited = [0 delimited length(confstr)+1];
    cpoi = 0;
    for k = 1:length(delimited)-1
        range = confstr(delimited(k)+1:delimited(k+1)-1);
        hyphen = strfind(range,'-');
        if length(hyphen) > 1
            exceptions = {MException('distance_distribution:wrong_conformer_range', 'Conformer range for site 2 contains more than one hyphen')};
            return
        end
        if isempty(hyphen) % single conformer
            cpoi = cpoi + 1;
            % process start/end syntax
            if strcmp(strtrim(range),'start')
                conformers2(cpoi) = 1;
            elseif strcmp(strtrim(range),'end')
                conformers2(cpoi) = length(entity.populations);
            else
                conformers2(cpoi) = str2double(range);
            end
        else % conformer range
            st_str = range(1:hyphen-1);
            % process start/end syntax
            if strcmp(strtrim(st_str),'start')
                st_val = 1;
            elseif strcmp(strtrim(st_str),'end')
                st_val = length(entity.populations);
            else
                st_val = str2double(st_str);
            end
            en_str = range(hyphen+1:end);
            % process start/end syntax
            if strcmp(strtrim(en_str),'start')
                en_val = 1;
            elseif strcmp(strtrim(en_str),'end')
                en_val = length(entity.populations);
            else
                en_val = str2double(en_str);
            end
            crange = st_val:en_val;
            conformers2(cpoi+1:cpoi+length(crange)) = crange;
            cpoi = cpoi + length(crange);
        end
    end
    conformers2 = conformers2(1:cpoi);
end

% for coupled mode, check whether the conformer ranges are consistent
if options.coupled
    if length(conformers1) ~= length(conformers2)
        exceptions = {MException('distance_distribution:inconsistent_conformer_ranges', 'Conformer ranges for sites have different length')};
        return
    end
    if max(abs(conformers1 - conformers2)) > 0
        exceptions = {MException('distance_distribution:inconsistent_conformer_ranges', 'Conformer selections for sites differ')};
        return
    end
end

% initialize empty distribution
distribution = zeros(size(r_axis));

profile on
% in coupled mode, loop over one conformer range
if options.coupled
    for kconf = conformers1
        site1 = sprintf('{%i}%s',kconf,residue1);
        [argsout,entity,exceptions] = get_label(entity,label1,{'positions','populations'},site1);
        if ~isempty(exceptions{1})
            distribution = [];
            return
        end
        positions1 = argsout{1}{1};
        populations1 = argsout{2}{1};
        site2 = sprintf('{%i}%s',kconf,residue2);
        [argsout,entity,exceptions] = get_label(entity,label2,{'positions','populations'},site2);
        if ~isempty(exceptions{1})
            distribution = [];
            return
        end
        positions2 = argsout{1}{1};
        populations2 = argsout{2}{1};
        % exclude, if one of the sites could not be labelled, raise an
        % exception
        if isempty(populations1) || isempty(populations2)
            warnings = warnings + 1;
            exceptions{warnings} = MException('distance_distribution:failed_site_pair', 'For %s-%s, the distance distribution is empty.',site1,site2);
            missing = 0; % because the warning was already raised
        else % if it worked out, add distribution
            [pair_distribution,missing] = get_pair_distribution(r_axis,...
                positions1,populations1,positions2,populations2,options.smoothing);
            if sum(isnan(pair_distribution))
                disp('Aber hallo!');
            end
            distribution = distribution + (1-missing)*pair_distribution;
        end
        if missing > 0.5e-2
            warnings = warnings + 1;
            exceptions{warnings} = MException('distance_distribution:missing_population', 'For %s-%s, %4.1f%% of population are outside distance range',site1,site2,100*missing);
        end
    end
% in uncoupled mode, loop over both conformer ranges
else
    % precompute the two sets of positions and populations
    all_positions1 = cell(1,length(conformers1));
    all_populations1 = cell(1,length(conformers1));
    for kconf1 = conformers1
        site1 = sprintf('{%i}%s',kconf1,residue1);
        [argsout,entity,exceptions] = get_label(entity,label1,{'positions','populations'},site1);
        if ~isempty(exceptions{1})
            distribution = [];
            return
        end
        all_positions1{kconf1} = argsout{1}{1};
        all_populations1{kconf1} = argsout{2}{1};
    end
    all_positions2 = cell(1,length(conformers2));
    all_populations2 = cell(1,length(conformers2));
    for kconf2 = conformers2
        site2 = sprintf('{%i}%s',kconf2,residue2);
        [argsout,entity,exceptions] = get_label(entity,label2,{'positions','populations'},site2);
        if ~isempty(exceptions{1})
            distribution = [];
            return
        end
        all_positions2{kconf2} = argsout{1}{1};
        all_populations2{kconf2} = argsout{2}{1};
    end
    % compute the sum of the distance distributions with pre-computed rotamers 
    for kconf1 = conformers1
        positions1 = all_positions1{kconf1};
        populations1 = all_populations1{kconf1};
        for kconf2 = conformers2
            positions2 = all_positions2{kconf2};
            populations2 = all_populations2{kconf2};
            if isempty(populations1) || isempty(populations2)
                warnings = warnings + 1;
                exceptions{warnings} = MException('distance_distribution:failed_site_pair', 'For %s-%s, the distance distribution is empty.',site1,site2);
            else
                [pair_distribution,missing] = get_pair_distribution(r_axis,...
                    positions1,populations1,positions2,populations2,options.smoothing);
                distribution = distribution + (1-missing)*pair_distribution;
                if missing > 0.5e-2
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('distance_distribution:missing_population', 'For %s-%s, %4.1f%% of population are outside distance range',site1,site2,100*missing);
                end
            end
        end
    end
end

% renormalize
distribution = distribution/sum(distribution);

% rescale units if requested
if strcmpi(options.units,'density')
    distribution = distribution/options.resolution;
end

if sum(isnan(distribution))
    distribution = [];
    exceptions = {MException('distance_distribution:empty_distribution', 'Distance distribution is empty for this site pair in the requested range')};
end

profile viewer

function [pair_distribution,missing] = get_pair_distribution(rax,positions1,populations1,positions2,populations2,smoothing)
% Computes a pair distance distribution as well as the fraction of missing
% population, the distribution is normalized to unity

%initialize output
pair_distribution = zeros(size(rax));

rmin = rax(1);
dr = rax(2)-rax(1);
nr = length(rax);

% normalize population vectors
populations1 = populations1/sum(populations1);
populations2 = populations2/sum(populations2);

[n_rot_1,~] = size(positions1);   % returns number of rotamers at position 1
[n_rot_2,~] = size(positions2);   % -//- at position 2

a2 = repmat(sum(positions1.^2,2),1,n_rot_2);
b2 = repmat(sum(positions2.^2,2),1,n_rot_1).';
pair_dist = sqrt(abs(a2 + b2 - 2*positions1*positions2.'));
weights = populations1*populations2.';
indices = reshape(pair_dist,1,n_rot_1*n_rot_2);
indices = 1 + round((indices-rmin)/dr);
weights = reshape(weights,1,n_rot_1*n_rot_2);
missing = sum(weights(indices > nr));
weights = weights(indices <= nr);
indices = indices(indices <= nr);
missing = missing + sum(weights(indices < 1));
weights = weights(indices >= 1);
indices = indices(indices >= 1);
for k = 1:length(indices)
    pair_distribution(indices(k)) = pair_distribution(indices(k)) + weights(k);
end

% Convolution with broadening function
vari = (0:dr:(nr-1)*dr)/(sqrt(2)*smoothing);
broadening = exp(-vari.^2); % broadening function
inverted_distribution = ifft(pair_distribution).*ifft(broadening);
pair_distribution = real(fft(inverted_distribution));

% normalization
pair_distribution = (1-missing)*pair_distribution/max(pair_distribution);

