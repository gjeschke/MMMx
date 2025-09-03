function [backbones,pop,exceptions] = get_backbones_ensemble(ensemble,chain,range,options)
%
% GET_BACKBONES_ENSEMBLE Retrieve backbone coordinates for an ensemble
%
%   [backbones,pop,exceptions] = GET_BACKBONES_ENSEMBLE(ensemble,chain,range)
%   Provides backbone coordinates and corresponding residue numbers for all
%   conformers of all biopolymer chains of an ensemble structure as well as
%   conformer populations and, if the function fails, exceptions,
%   optionally, retrieval can be confined to a chain and a residue range
%   within the chain
%
% INPUT
% ensemble  either an entity or a file name specifying the ensemble with
%           possible extensions
%           .ens   ensemble definition file with populations
%           .pdb   PDB file, filename can contain wildcards, all files 
%                  matching it are processed, 
%           none   extension .pdb is appended
% chain     chain identifier, defaults to all chains, can be a string for
%           several chains
% range     cell vectôr of residue ranges, either double(1,2), or vector of
%           negative residue numbers,
%           if empty or missing, it defaults to all residues
% options   .full   flag, if true, N and C coordinates for amino acid
%                   residues are also extracted
%
% OUTPUT
% backbones     struct with fields (chain) for all chains, each chain has
%               fields, empty array if function fails
%               not a biopolymer
%               .type 1 for peptide, 2 for nucleic acid, for mixed chains,
%                     the type with more backbone atoms prevails
%               .mono (1,nm) double, numbers of the residues
%               .bb   (1,c) cell of (nm,3) coordinate array for CA (peptide)
%                     or C4' (nucleic acid) atoms
%               .slc  sequence in single-letter code
%               if options.full is true
%               .N    (1,c) cell of (nm,3) coordinate array for N atoms
%               .C    (1,c) cell of (nm,3) coordinate array for C atoms
% pop           (1,C) population vector for C conformers
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
%
% exceptions occur if chain length or residue identity is inconsistent
% between conformers in any chain
% however, if all biopolymers in the first PDB file are encountered in
% later files and are consistent, additional biopolymer chains in later
% files do not raise an exception

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
backbones = [];
exceptions = {[]};
pop = [];
popdef = [];

% set defaults for missing input
if ~exist('chain','var') || isempty(chain)
    all_chains = true;
    selected_chain = '';
else
    all_chains = false;
    selected_chain = chain;
end
if ~exist('range','var') || isempty(range)
    if isempty(selected_chain)
        ranges = cell(1,1);
        ranges{1} = [1, 1e10];
    else
        ranges = cell(length(selected_chain),1);
        for k = 1:length(selected_chain)
            ranges{k} = [1, 1e10];
        end
    end
else
    ranges = range;
end

if ~exist('options','var') || isempty(options) || ~isfield(options,'full')
    options.full = false;
end

if isstruct(ensemble) % input is an entity
    pop = ensemble.populations;
    for c = 1:length(ensemble.populations) % loop over all conformers
        if c == 1 % initialize chains
            chains = fieldnames(ensemble);
            for kc = 1:length(chains)
                chain = chains{kc};
                if ~all_chains 
                    index = strfind(selected_chain,chain);
                    if isempty(index)
                        continue
                    end
                    range = ranges{index};
                end
                if isstrprop(chain(1),'upper') % chain fields start with a capital
                    backbone = get_backbone(ensemble,chain,1,range,options);
                    if isempty(backbone.aa) && isempty(backbone.nt) % skip chains that are not a biopolymer
                        continue
                    end
                    backbones.(chain).bb = cell(1,length(ensemble.populations));
                    if length(backbone.aa) >= length(backbone.nt) % peptide prevails
                        backbones.(chain).type = 1;
                        backbones.(chain).mono = backbone.aa;
                        backbones.(chain).bb{1} = backbone.CA;
                        backbones.(chain).slc = backbone.res;
                        if options.full
                            backbones.(chain).N{1} = backbone.N;
                            backbones.(chain).C{1} = backbone.C;
                        end
                    else % its a nucleic acid
                        backbones.(chain).type = 2;
                        backbones.(chain).mono = backbone.nt;
                        backbones.(chain).bb{1} = backbone.C4p;
                        backbones.(chain).slc = backbone.nuc;
                    end
                end
            end
        else % not the first conformer
            for kc = 1:length(chains)
                chain = chains{kc};
                if ~all_chains 
                    index = strfind(selected_chain,chain);
                    if isempty(index)
                        continue
                    end
                    range = ranges{index};
                end
                if isstrprop(chain(1),'upper') % chain fields start with a capital
                    backbone = get_backbone(ensemble,chain,c,range,options);
                    switch backbones.(chain).type
                        case 1 % peptide
                            if length(backbones.(chain).mono) ~= length(backbone.aa)
                                exceptions = {MException('get_backbones_ensemble:chain_length_mismatch',...
                                    'Conformer %i of chain %s differs in length from first conformer',c,chain)};
                                backbones = [];
                                return
                            end
                            deviation = sum(abs(backbones.(chain).mono - backbone.aa));
                            if deviation > 0
                                exceptions = {MException('get_backbones_ensemble:chain_index_mismatch',...
                                    'Conformer %i of chain %s differs in residue numbers from first conformer',c,chain)};
                                backbones = [];
                                return
                            end
                            backbones.(chain).bb{c} = backbone.CA;
                            if options.full
                                backbones.(chain).N{c} = backbone.N;
                                backbones.(chain).C{c} = backbone.C;
                            end
                        case 2 % nucleic acid
                            if length(backbones.(chain).mono) ~= length(backbone.nt)
                                exceptions = {MException('get_backbones_ensemble:chain_length_mismatch',...
                                    'Conformer %i of chain %s differs in length from first conformer',c,chain)};
                                backbones = [];
                                return
                            end
                            deviation = sum(abs(backbones.(chain).mono - backbone.nt));
                            if deviation > 0
                                exceptions = {MException('get_backbones_ensemble:chain_index_mismatch',...
                                    'Conformer %i of chain %s differs in residue numbers from first conformer',c,chain)};
                                backbones = [];
                                return
                            end
                            backbones.(chain).bb{c} = backbone.C4p;
                    end
                end
            end
        end
    end
else % file name is given
    [~,~,ext] = fileparts(ensemble);
    if isempty(ext)
        ensemble = strcat(ensemble,'.pdb');
        all_files = dir(ensemble); % find all files that match the pattern
    elseif strcmp(ext,'.ens')
        [all_files,popdef] = rd_ensemble_definition(ensemble);
    elseif ~strcmp(ext,'.pdb')
        exceptions = {MException('get_backbones_ensemble:undefined_input_file',...
            'Cannot determine file type for extension %s',ext)};
        return
    else
        all_files = dir(ensemble); % find all files that match the pattern
    end
    for f = 1:length(all_files)
        entity = get_pdb(all_files(f).name);
        if f == 1 % initialize output
            % we can do this by a recursive algorithm
            [backbones,pop,exceptions] = get_backbones_ensemble(entity,selected_chain,range);
            if ~isempty(exceptions{1})
                return
            end
        else % these are additional conformers from more PDB files
            [backbones2,pop2,exceptions] = get_backbones_ensemble(entity,selected_chain,range);
            if ~isempty(exceptions{1})
                exceptions{2} = MException('get_backbones_ensemble:wrong_file',...
                    'File %s does not match expected format',all_files(f).name);
                backbones = [];
                return
            end
            chains = fieldnames(backbones);
            pop = [pop pop2]; %#ok<AGROW>
            for kc = 1:length(chains) % try all biopolymer chains present in PDB file 1
                chain = chains{kc};
                if ~isfield(backbones2,chain) % this file does not contain a biopolymer that occurred in the first file
                    exceptions = {MException('get_backbones_ensemble:chain_missing',...
                        'File %s misses chain %s',all_files(f).name,chain)};
                    backbones = [];
                    return
                end
                % check for consistency of this chain with the one in the first PDB file
                if backbones.(chain).type ~= backbones2.(chain).type
                    exceptions = {MException('get_backbones_ensemble:type_mismatch',...
                        'In file %s chain %s has type %i, but should have type %i',...
                        all_files(f).name,chain,backbones2.(chain).type,backbones.(chain).type)};
                    backbones = [];
                    return
                end
                if length(backbones.(chain).mono) ~= length(backbones2.(chain).mono)
                    exceptions = {MException('get_backbones_ensemble:length_mismatch',...
                        'In file %s chain %s has length %i, but should have length %i',...
                        all_files(f).name,chain,length(backbones2.(chain).mono),length(backbones.(chain).mono))};
                    backbones = [];
                    return
                end
                known_conformers = length(backbones.(chain).bb); % how many conformers have already been read?
                % add the new conformers
                for c = 1:length(backbones2.(chain).bb)
                    backbones.(chain).bb{known_conformers+c} = backbones2.(chain).bb{c};
                    if options.full
                        backbones.(chain).N{known_conformers+c} = backbones2.(chain).N{c};
                        backbones.(chain).C{known_conformers+c} = backbones2.(chain).C{c};
                    end
                end
            end
        end
    end
end

if ~isempty(popdef)
    pop = popdef;
end
% renormalize populations
pop = pop/sum(pop);
