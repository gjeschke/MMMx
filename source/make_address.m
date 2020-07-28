function [address,exceptions] = make_address(chain,residue,atom,conformer,location)
%
% MAKE_ADDRESS Make MMMx address from components
%
%   [address,exceptions] = MAKE_ADDRESS
%   Combines object information into an MMMx address
%
%   [address,exceptions] = MAKE_ADDRESS(chain,residue,atom,conformer,location)
%   Analyzes object information and combines it into an MMMx address
%
% INPUT
% chain     chain identifier, string, preferably single capital letter or
%           cell of strings
% residue   vector of integer residue numbers or a string (residue name)
%           or a cell string (residue names)
% atom      atom names, string or cell of strings, each string preferably 
%           not longer than four characters
% conformer integer or vector of integers
% location  atom location tag, string or cell of strings, each string 
%           preferably a single capital letter
%
% OUTPUT
% address      MMMx address string, empty, if an error occurred
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

address = '';
exceptions{100} = [];
warnings = 0; % counter for warnings

if exist('chain','var') && ~isempty(chain)
    % check chain argument for valid type
    if iscell(chain)
        if ~ischar(chain{1})
            exceptions = {MException('make_address:chainID_no_char', 'Chain identifier is not a character or string')};
            return
        else
            chains = chain{1};
            for k = 2:length(chain)
                if ischar(chain{k})
                    chains = strcat(chains,sprintf(',%s',chain{k}));
                else
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('make_address:part_of_chainID_no_char', 'Chain identifier is not completely a character or string');
                end
            end
        end
    elseif ~ischar(chain)
        exceptions = {MException('make_address:chainID_no_char', 'Chain identifier is not a character or string')};
        return
    else
        chains = chain;
    end
    
    if ~isempty(chains)
        address = sprintf('(%s)',chains);
    end
end

conformers = '';
% translate vector of conformer numbers to abbreviated list of ranges
if exist('conformer','var') && ~isempty(conformer)
    if ischar(conformer)
        conformers = conformer;
    else
        conformers = sprintf('%i',conformer(1));
        cres = conformer(1);
        continuous = false;
        for k = 2:length(conformer)
            if conformer(k) ~= cres + 1
                if continuous
                    conformers = sprintf('%s-%i',conformers,cres);
                end
                conformers = sprintf('%s,%i',conformers,conformer(k));
                continuous = false;
            else
                continuous = true;
            end
            cres = conformer(k);
        end
        if continuous
            conformers = sprintf('%s-%i',conformers,conformer(k));
        end
    end
end

if ~isempty(conformers)
    address = sprintf('%s{%s}',address,conformers);
end

residues = '';
% translate vector of residue numbers to abbreviated list of ranges
if exist('residue','var') && ~isempty(residue)
    if iscell(residue)
        if ~ischar(residue{1})
            exceptions = {MException('make_address:residue_name_no_char', 'First residue name in list is not a string')};
            address = '';
            return
        else
            residues = upper(residue{1});
            for k = 2:length(residue)
                if ischar(residue{k})
                    residues = strcat(residues,sprintf(',%s',upper(residue{k})));
                else
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('make_address:part_of_residue_name_list_no_char', 'Residue name list contains non-strings');
                end
            end
        end
    elseif ischar(residue)
        residues = upper(residue);
    else
        residues = sprintf('%i',residue(1));
        cres = residue(1);
        continuous = false;
        for k = 2:length(residue)
            if residue(k) ~= cres + 1
                if continuous
                    residues = sprintf('%s-%i',residues,cres);
                end
                residues = sprintf('%s,%i',residues,residue(k));
                continuous = false;
            else
                continuous = true;
            end
            cres = residue(k);
        end
        if continuous
            residues = sprintf('%s-%i',residues,residue(k));
        end
    end
end

if ~isempty(residues)
    address = strcat(address,residues);
end

atoms = '';
if exist('atom','var') && ~isempty(atom)
    if iscell(atom)
        if ~ischar(atom{1})
            exceptions = {MException('make_address:atomID_no_char', 'atom identifier is not a character or string')};
        else
            atoms = atom{1};
            for k = 2:length(atom)
                if ischar(atom{k})
                    atoms = strcat(atoms,sprintf(',%s',atom{k}));
                else
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('make_address:part_of_atomID_no_char', 'atom identifier is not completely a character or string');
                end
            end
        end
    elseif ~ischar(atom)
        exceptions = {MException('make_address:atomID_no_char', 'atom identifier is not a character or string')};
        return
    else
        atoms = atom;
    end
end

if ~isempty(atoms)
    address = sprintf('%s.%s',address,atoms);
end

locations = '';
if exist('location','var') && ~isempty(location) 
    if iscell(location)
        if ~ischar(location{1})
            exceptions = {MException('make_address:locationID_no_char', 'location identifier is not a character or string')};
        else
            locations = location{1};
            for k = 2:length(location)
                if ischar(location{k})
                    locations = strcat(locations,sprintf(',%s',location{k}));
                else
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('make_address:part_of_locationID_no_char', 'location identifier is not completely a character or string');
                end
            end
        end
    elseif ~ischar(location)
        exceptions = {MException('make_address:locationID_no_char', 'location identifier is not a character or string')};
        return
    else
        locations = location;
    end
end

if ~isempty(locations)
    address = sprintf('%s:%s',address,locations);
end

exceptions = exceptions(1:warnings);