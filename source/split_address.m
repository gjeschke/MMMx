function [chains,residues,atoms,conformers] = split_address(address)
%
% SPLIT_ADDRESS Split MMMx address into ranges of objects
%
%   [chains,residues,atoms,conformers] = split_address(address)
%   Derives untested object ranges from an MMMx address 
%
% INPUT:
% 
% address   MMMx address string, such as (A,C,E-G){3,7}182-185.CA,N
%
% OUTPUT:
%
% chains        cell (1,c) of chain tags, example: {'A','C','E','F','G'}
% residues      double(1,r) of residue numbers, example [182,183,184,185]
% atoms         cell(1,a) of atom identifiers, example {'CA','N'} 
% conformers    double(1,C) of conformer numbers, example [3,7]
%
% range identifier '-' is allowed only on chain, residue, and conformer
% level
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% initialize empty output

chains = '';
residues = [];
atoms = {};
conformers = [];

% split address into parts
pa = strfind(address,'(');
pe = strfind(address,')');
if ~isempty(pa) && ~isempty(pe)
    chain_string = address(pa+1:pe-1);
    address = [address(1:pa-1) address(pe+1:end)];
else
    chain_string = '';
end

pa = strfind(address,'{');
pe = strfind(address,'}');
if ~isempty(pa) && ~isempty(pe)
    conformer_string = address(pa+1:pe-1);
    address = [address(1:pa-1) address(pe+1:end)];
else
    conformer_string = '';
end

parts = split(address,'.');
if ~isempty(parts{1})
    residue_string = parts{1};
else
    residue_string = '';
end
if length(parts) > 1
    atom_string = parts{2};
else
    atom_string = '';
end

% analyze substrings

if ~isempty(chain_string)
    chain_args = split(chain_string,',');
end
for k = 1:length(chain_args)
    if contains(chain_args{k},'-')
        range = split(chain_args{k},'-');
        chains = strcat(chains,char([double(range{1}):double(range{2})])); %#ok<NBRAK>
    else
        chains = strcat(chains,chain_args{k});
    end
end

if ~isempty(residue_string)
    residue_args = split(residue_string,',');
else
    residue_args = {};
end
for k = 1:length(residue_args)
    if contains(residue_args{k},'-')
        range = split(residue_args{k},'-');
        residues = [residues,str2double(range{1}):str2double(range{2})]; %#ok<AGROW>
    else
        residues = [residues,str2double(residue_args{k})]; %#ok<AGROW>
    end
end

if ~isempty(conformer_string)
    conformer_args = split(conformer_string,',');
else
    conformer_args = {};
end
for k = 1:length(conformer_args)
    if contains(conformer_args{k},'-')
        range = split(conformer_args{k},'-');
        conformers = [conformers,str2double(range{1}):str2double(range{2})]; %#ok<AGROW>
    else
        conformers = [conformers,str2double(conformer_args{k})]; %#ok<AGROW>
    end
end

if ~isempty(atom_string)
    atoms = split(atom_string,',');
end
