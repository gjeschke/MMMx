function [r_axes,trivariate,entity,exceptions] = trivariate_distribution(entity,site1,label1,site2,label2,site3,label3,options)
%
% TRIVARIATE_DISTRIBUTION Trivariate distance distribution between labels 
%                         or atoms
%
%   [r_axes,trivariate,entity,exceptions] =
%   TRIVARIATE_DISTRIBUTION(entity,site1,label1,site2,label2,site3,label3)
%   Trivariate distance distribution between residues site1, site2, and 
%   site 3 with label/atom specifications label1, label2, and label3
%
%   [entity,exceptions] =
%   TRIVARIATE_DISTRIBUTION(entity,site1,label1,site2,label2,site3,label3,options)
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
% site3     residue address string for site 3 with optional conformer 
%           specifier, use '{*}' for all conformers, defaults to '{1}'
% label3    label name, such as 'mtsl' or 'I1A' or atom specifier for this
%           residue, such as 'atom.CA', for label at site 3
% options   computation options
%           .rmin         minimum distance (Angstroem), default 10, double
%           .rmax         maximum distance (Angstroem), default 150, double
%           .resolution   resolution (Angstroem), default 0.5, double
%           .units        string, 'probability' makes sum(distribution) = 1
%                         'density' returns probability density in units of
%                         (1/Angstroem), default: 'density'
%
% OUTPUT
% r_axes       cell(1,3) of distance axes in units of Angstroem, sequence
%              is site1-site2, site1-site3, site2-site3
% trivariate   trivariate distance distribution
% entity       the input entity augmented by newly computed labels
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2022: Gunnar Jeschke

% initialize output
trivariate = [];

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

r_axis = options.rmin:options.resolution:options.rmax;
r_axes{1} = r_axis;
r_axes{2} = r_axis;
r_axes{3} = r_axis;

% extract conformer selection of site1 if any, set default, if none
si = strfind(site1,'{');
ei = strfind(site1,'}');
if ~isempty(si) && ~isempty(ei)
    confstr = site1(si+1:ei-1);
    residue1 = [site1(1:si-1) site1(ei+1:end)];
else
    confstr = '*';
    residue1 = site1;
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
            exceptions = {MException('trivariate_distribution:wrong_conformer_range', 'Conformer range for site 1 contains more than one hyphen')};
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

% extract conformer selection of site2 if any, set default, if none
si = strfind(site2,'{');
ei = strfind(site2,'}');
if ~isempty(si) && ~isempty(ei)
    confstr = site2(si+1:ei-1);
    residue2 = [site2(1:si-1) site2(ei+1:end)];
else
    confstr = '*';
    residue2 = site2;
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
            exceptions = {MException('trivariate_distribution:wrong_conformer_range', 'Conformer range for site 2 contains more than one hyphen')};
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

% extract conformer selection of site3 if any, set default, if none
si = strfind(site3,'{');
ei = strfind(site3,'}');
if ~isempty(si) && ~isempty(ei)
    confstr = site3(si+1:ei-1);
    residue3 = [site3(1:si-1) site3(ei+1:end)];
else
    confstr = '*';
    residue3 = site3;
end
if strcmp(strtrim(confstr),'*')
    conformers3 =  1:length(entity.populations);
    confstr = '';
end
% make conformer vector
if ~isempty(confstr)
    conformers3 = zeros(1,66000);
    delimited = strfind(confstr,',');
    delimited = [0 delimited length(confstr)+1];
    cpoi = 0;
    for k = 1:length(delimited)-1
        range = confstr(delimited(k)+1:delimited(k+1)-1);
        hyphen = strfind(range,'-');
        if length(hyphen) > 1
            exceptions = {MException('trivariate_distribution:wrong_conformer_range', 'Conformer range for site 2 contains more than one hyphen')};
            return
        end
        if isempty(hyphen) % single conformer
            cpoi = cpoi + 1;
            % process start/end syntax
            if strcmp(strtrim(range),'start')
                conformers3(cpoi) = 1;
            elseif strcmp(strtrim(range),'end')
                conformers3(cpoi) = length(entity.populations);
            else
                conformers3(cpoi) = str2double(range);
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
            conformers3(cpoi+1:cpoi+length(crange)) = crange;
            cpoi = cpoi + length(crange);
        end
    end
    conformers3 = conformers3(1:cpoi);
end

% for coupled mode, check whether the conformer ranges are consistent
if options.coupled
    if length(conformers1) ~= length(conformers2) || length(conformers1) ~= length(conformers3)
        exceptions = {MException('distance_distribution:inconsistent_conformer_ranges', 'Conformer ranges for sites have different length')};
        return
    end
    if max(abs(conformers1 - conformers2)) > 0 || max(abs(conformers1 - conformers3)) > 0
        exceptions = {MException('distance_distribution:inconsistent_conformer_ranges', 'Conformer selections for sites differ')};
        return
    end
end

% initialize empty distribution
trivariate = zeros(length(r_axes{1}),length(r_axes{2}),length(r_axes{3}));

% in coupled mode, loop over one conformer range
if options.coupled
    for kconf = conformers1
        site1 = sprintf('{%i}%s',kconf,residue1);
        [argsout,entity,exceptions] = get_label(entity,label1,{'positions','populations'},site1);
        if ~isempty(exceptions{1})
            trivariate = [];
            return
        end
        positions1 = argsout{1}{1};
        populations1 = argsout{2}{1};
        site2 = sprintf('{%i}%s',kconf,residue2);
        [argsout,entity,exceptions] = get_label(entity,label2,{'positions','populations'},site2);
        if ~isempty(exceptions{1})
            trivariate = [];
            return
        end
        positions2 = argsout{1}{1};
        populations2 = argsout{2}{1};
        site3 = sprintf('{%i}%s',kconf,residue3);
        [argsout,entity,exceptions] = get_label(entity,label3,{'positions','populations'},site3);
        if ~isempty(exceptions{1})
            trivariate = [];
            return
        end
        positions3 = argsout{1}{1};
        populations3 = argsout{2}{1};
        % exclude, if one of the sites could not be labelled, raise an
        % exception
        if isempty(populations1) || isempty(populations2) || isempty(populations3)
            warnings = warnings + 1;
            exceptions{warnings} = MException('trivariate_distribution:failed_site_triple', 'For %s-%s-%s, the trivariate distribution is empty.',site1,site2,site3);
            missing = 0; % because the warning was already raised
        else % if it worked out, add distribution
            [trivariate_c,missing] = get_trivariate_distribution(r_axes,...
                positions1,populations1,positions2,populations2,...
                positions3,populations3);
            if sum(isnan(trivariate_c))
                warnings = warnings + 1;
                exceptions{warnings} = MException('trivariate_distribution:failed_site_triple', 'For %s-%s-%s, the trivariate distribution contains NaN.',site1,site2,site3);
            else
                trivariate = trivariate + (1-missing)*trivariate_c;
            end
        end
        if missing > 0.5e-2
            warnings = warnings + 1;
            exceptions{warnings} = MException('trivariate_distribution:missing_population', 'For %s-%s-%s, %4.1f%% of population are outside distance range',site1,site2,site3,100*missing);
        end
    end
% no uncoupled mode implemented for trivariate distributions
else
    exceptions = {MException('trivariate_distribution:no_uncoupled_mode', 'For trivariate distributions, uncoupled mode is not implemented')};
    return
end

% renormalize
trivariate = trivariate/sum(sum(sum(trivariate)));

% rescale units if requested
if strcmpi(options.units,'density')
    trivariate = trivariate/options.resolution^3;
end

if sum(isnan(trivariate))
    trivariate = [];
    exceptions = {MException('trivariate_distribution:empty_distribution', 'Trivariate distribution is empty for this site triple in the requested range')};
end


function [trivariate,missing] = get_trivariate_distribution(r_axes,positions1,populations1,positions2,populations2,positions3,populations3)
% Computes a pair distance distribution as well as the fraction of missing
% population, the distribution is normalized to unity

%initialize output
trivariate = zeros(length(r_axes{1}),length(r_axes{2}),length(r_axes{3}));


rmin1 = r_axes{1}(1);
dr1 = r_axes{1}(2)-r_axes{1}(1);
nr1 = length(r_axes{1});
rmin2 = r_axes{2}(1);
dr2 = r_axes{2}(2)-r_axes{2}(1);
nr2 = length(r_axes{2});
rmin3 = r_axes{3}(1);
dr3 = r_axes{3}(2)-r_axes{3}(1);
nr3 = length(r_axes{3});

% normalize population vectors
populations1 = populations1/sum(populations1);
populations2 = populations2/sum(populations2);
populations3 = populations3/sum(populations3);

[n_rot_1,~] = size(positions1);   % returns number of rotamers at position 1
[n_rot_2,~] = size(positions2);   % -//- at position 2
[n_rot_3,~] = size(positions3);   % -//- at position 3

missing = 0;
for k1 = 1:n_rot_1
    weight0 = populations1(k1);
    for k2 = 1:n_rot_2
        weight1 = weight0 * populations2(k2);
        a = norm(positions1(k1,:)-positions2(k2,:));
        index_a = 1 + round((a-rmin1)/dr1);
        if index_a < 1 || index_a > nr1
            missing = missing + weight1;
        else
            for k3 = 1:n_rot_3
                weight = weight1 * populations3(k3);
                b = norm(positions1(k1,:)-positions3(k3,:));
                c = norm(positions2(k2,:)-positions3(k3,:));
                index_b = 1 + round((b-rmin2)/dr2);
                index_c = 1 + round((c-rmin3)/dr3);
                if index_b < 1 || index_b > nr2 || index_c < 1 || index_c > nr3
                    missing = missing + weight;
                else
                    trivariate(index_a,index_b,index_c) = trivariate(index_a,index_b,index_c) + weight;
                end
            end
        end
    end
end

% smoothing
trivariate = smooth3(trivariate,'gaussian',5,2);
% normalization
trivariate = (1-missing)*trivariate/max(max(max(trivariate)));


