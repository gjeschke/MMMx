function [addresses,exceptions] = cx_to_mmmx_address(spec)
%
% CX_TO_MMMX_ADDRESS Translate ChimeraX target specification
%
%   [addresses,exceptions] = CX_TO_MMMX_address(spec)
%   Translates a ChimeraX target specification into a set of MMMx
%   addresses and ChimeraX structure identifiers, only a subset of ChimeraX
%   syntax is supported, use cx_select, if you need more functionality
%
% INPUT
% spec      ChimeraX target specification, see:
%           https://www.cgl.ucsf.edu/chimerax/docs/user/commands/atomspec.html
%
% OUTPUT
% addresses  array of structures with fields .id and .address
%            .id        structure identifier (integer) in ChimeraX
%                       empty, if not specified
%            .address   MMMx address string, empty if the ChimeraX
%                       specifier could not be interpreted
%            .conformer conformer number
%            .chain     chain identifier
%            .residue   residue identifier
%            .atom      atom identifier
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initalize empty output
addresses(1).id = [];
addresses(1).address = '';
exceptions = {[]};

id_done = false;
chain_done = false;
residue_done = false;
atom_str = '';

% everything before a hash is discarded
poi = strfind(spec,'#');
if ~isempty(poi)
    remainder = spec(poi+1:end);
else
    id_done = true;
    id_str = '';
    remainder = spec;
end

% is a chain specification present?
poi = strfind(remainder,'/');
if ~isempty(poi)
    if ~id_done
        id_str = remainder(1:poi-1);
        id_done = true;
    end
    remainder = remainder(poi+1:end);
else
    chain_done = true;
    chain_str = '';
end

% is a residue specification present?
poi = strfind(remainder,':');
if ~isempty(poi)
    if ~id_done
        id_str = remainder(1:poi-1);
        id_done = true;
    end
    if ~chain_done
        chain_str = remainder(1:poi-1);
        chain_done = true;
    end
    remainder = remainder(poi+1:end);
else
    residue_done = true;
    residue_str = '';
end

% is an atom specification present?
poi = strfind(remainder,'@');
if ~isempty(poi)
    if ~id_done
        id_str = remainder(1:poi-1);
        id_done = true;
    end
    if ~chain_done
        chain_str = remainder(1:poi-1);
        chain_done = true;
    end
    if ~residue_done
        residue_str = remainder(1:poi-1);
        residue_done = true;
    end
    atom_str = remainder(poi+1:end);
    remainder = '';
end

if ~id_done
    id_str = remainder;
    remainder = '';
end
if ~chain_done
    chain_str = remainder;
    remainder = '';
end
if ~residue_done
    residue_str = remainder;
    remainder = '';
end

if ~isempty(remainder)
    exceptions = {MException('cx_to_mmmx_address:surplus_specifier_part', 'Surplus part of specifier: %s', remainder)};
    return
end

% check for extended ChimeraX syntax that is not supported
valid = check_part(id_str,'.');
if ~valid
    exceptions = {MException('cx_to_mmmx_address:invalid_id_syntax', 'Invalid syntax of ID: %s', id_str)};
    return
end

% separate ID string, if a list
ids = {[]};
conformers = {''};
if ~isempty(id_str)
    delimited = strfind(id_str,',');
    ids = cell(1,length(delimited)+1);
    conformers = cell(1,length(delimited)+1);
    delimited = [0 delimited length(id_str)+1];
    for k = 1:length(delimited)-1
        id_sub_str = id_str(delimited(k)+1:delimited(k+1)-1);
        % separate identifier string into ID and conformers, if any specified
        poi = strfind(id_sub_str,'.');
        if ~isempty(poi)
            conformers{k} = id_sub_str(poi+1:end);
            ids{k} = str2double(id_sub_str(1:poi-1));
        else
            conformers{k} = '';
            ids{k} = str2double(id_sub_str);
        end
    end
end

% check for extended ChimeraX syntax that is not supported
for k = 1:length(conformers)
    conf_str = conformers{k};
    valid = check_part(conf_str);
    if ~valid
        exceptions = {MException('cx_to_mmmx_address:invalid_conf_syntax', 'Invalid syntax of conformer part: %s', conf_str)};
        return
    end
end
valid = check_part(chain_str);
if ~valid
    exceptions = {MException('cx_to_mmmx_address:invalid_chain_syntax', 'Invalid syntax of chain part: %s', chain_str)};
    return
end
valid = check_part(residue_str);
if ~valid
    exceptions = {MException('cx_to_mmmx_address:invalid_residue_syntax', 'Invalid syntax of residue part: %s', residue_str)};
    return
end
valid = check_part(atom_str);
if ~valid
    exceptions = {MException('cx_to_mmmx_address:invalid_atom_syntax', 'Invalid syntax of atom part: %s', atom_str)};
    return
end

addresses(length(ids)).id = [];
addresses(length(ids)).address = '';
for k = 1:length(ids)
    addresses(k).id = ids{k};
    address = '';
    conf_str = conformers{k};
    if ~isempty(conf_str)
        address = sprintf('{%s}',conf_str);
    end
    if ~isempty(chain_str)
        address = sprintf('%s(%s)',address,chain_str);
    end
    if ~isempty(residue_str)
        address = sprintf('%s%s',address,residue_str);
    end
    if ~isempty(atom_str)
        address = sprintf('%s.%s',address,atom_str);
    end
    addresses(k).address = address;
    addresses(k).conformer = conformers{k};
    addresses(k).chain = chain_str;
    addresses(k).residue = residue_str;
    addresses(k).atom = atom_str;
end
% fprintf(1,'chain part    : |%s|\n',chain_str);
% fprintf(1,'residue part  : |%s|\n',residue_str);
% fprintf(1,'atom part     : |%s|\n',atom_str);

function valid = check_part(part,allowed)
% checks whether a string is valid for translation

valid = true;
if contains(part,'#') && ~contains(allowed,'#')
    valid = false;
    return
end
if contains(part,'/') && ~contains(allowed,'/')
    valid = false;
    return
end
if contains(part,':') && ~contains(allowed,':')
    valid = false;
    return
end
if contains(part,'@') && ~contains(allowed,'@')
    valid = false;
    return
end
if contains(part,'?') && ~contains(allowed,'?')
    valid = false;
    return
end
if contains(part,'|') && ~contains(allowed,'|')
    valid = false;
end
if contains(part,'=') && ~contains(allowed,'=')
    valid = false;
    return
end
if contains(part,'>') && ~contains(allowed,'>')
    valid = false;
    return
end
if contains(part,'<') && ~contains(allowed,'<')
    valid = false;
    return
end
if contains(part,'~') && ~contains(allowed,'~')
    valid = false;
    return
end
if contains(part,'+') && ~contains(allowed,'+')
    valid = false;
    return
end
if contains(part,'[') && ~contains(allowed,'[')
    valid = false;
    return
end
if contains(part,']') && ~contains(allowed,']')
    valid = false;
    return
end
if contains(part,'.') && ~contains(allowed,'.') % to avoid more than two levels
    valid = false;
    return
end