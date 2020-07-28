function [entity,exceptions] = select(entity,address,overwrite,unselect)
%
% SELECT Select objects in an entity (MMMx internal) 
%
%   [entity,exceptions] = SELECT
%   Selects all objects in an entity that match an MMMx address
%
%   [entity,exceptions] = SELECT(entity,address,overwrite,unselect)
%   Expands an MMMx address into vectors or field name lists of matching 
%   objects at all levels and marks them as selected in the entity 
%
% INPUT
% entity    entity in an MMMx format
% address   MMMx address string, if a location is given and atoms with
%           empty location tag should be selected, this must be explicit in
%           the address
% overwrite optional Boolean flag, if true, an existing selection is
%           deleted, if false, the new selection adds to an existing one,
%           defaults to false, overwrite with empty address unselects all
% unselect  optional Boolean flag, if true, the specified objects are
%           unselected, defaults to false, true unselect forces overwrite
%           to be false
%
% OUTPUT
% entity       the input entity with the selections made (or changed)
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

profile on

exceptions{100} = [];
warnings = 0; % counter for warnings

% initialize overwrite, if it does not exist
if ~exist('overwrite','var') || isempty(overwrite)
    overwrite = false;
end

% initialize invert, if it does not exist
if ~exist('unselect','var') || isempty(unselect)
    unselect = false;
elseif invert
    overwrite = false;
end

if overwrite % delete existing selection
    entity.selected = [];
    chains = fieldnames(entity);
    for kc = 1:length(chains)
        chain = chains{kc};
        if isstrprop(chain(1),'upper') % chain fields start with a capital
            entity.(chain).selected = 0;
            residues = fieldnames(entity.(chain));
            for kr = 1:length(residues)
                residue = residues{kr};
                if strcmp(residue(1),'R') % these are residue fields
                    entity.(chain).(residue).selected = [];
                    atoms = fieldnames(entity.(chain).(residue));
                    for ka = 1:length(atoms)
                        atom = atoms{ka};
                        if isstrprop(atom(1),'upper') % these are atom fields
                            entity.(chain).(residue).(atom).selected = [];
                        end
                    end
                end
            end
        end
    end
end

% expand address
% separate chain part and divide into single chain tags
chains = '';
si = strfind(address,'(');
ei = strfind(address,')');
if ~isempty(si) && ~isempty(ei)
    chains = address(si+1:ei-1);
end
delimited = strfind(chains,',');
chain_fields = cell(1,length(delimited)+1);
delimited = [0 delimited length(chains)+1];
for k = 1:length(chain_fields)
    chain_fields{k} = chains(delimited(k)+1:delimited(k+1)-1);
end
% remove part of address that was already processed
if ~isempty(ei)
    address = address(ei+1:end);
end

% separate conformer part and make conformer vector
conformers = zeros(1,66000);
si = strfind(address,'{');
ei = strfind(address,'}');
if ~isempty(si) && ~isempty(ei)
    confstr = address(si+1:ei-1);
end
delimited = strfind(confstr,',');
delimited = [0 delimited length(confstr)+1];
cpoi = 0;
for k = 1:length(delimited)-1
    range = confstr(delimited(k)+1:delimited(k+1)-1);
    hyphen = strfind(range,'-');
    if length(hyphen) > 1
        exceptions{1} = MException('select:wrong_conformer_range', 'Conformer range contains more than one hyphen');
        return
    end
    if isempty(hyphen) % single conformer
        cpoi = cpoi + 1;
        conformers(cpoi) = str2double(range);
    else % conformer range
        crange = str2double(range(1:hyphen-1)):str2double(range(hyphen+1:end));
        conformers(cpoi+1:cpoi+length(crange)) = crange;
        cpoi = cpoi + length(crange);
    end
end
conformers = conformers(1:cpoi);
% remove part of address that was already processed
if ~isempty(ei)
    address = address(ei+1:end);
end

% separate residues part and make residue vector and residue type
% cellstring
residues = zeros(1,66000);
rpoi = 0;
residue_types = cell(1,100);
rtpoi = 0;
ei = strfind(address,'.');
if ~isempty(ei)
    resstr = address(1:ei-1);
end
delimited = strfind(resstr,',');
delimited = [0 delimited length(resstr)+1];
for k = 1:length(delimited)-1
    range = resstr(delimited(k)+1:delimited(k+1)-1);
    hyphen = strfind(range,'-');
    if length(hyphen) > 1
        exceptions{1} = MException('select:wrong_conformer_range', 'Conformer range contains more than one hyphen');
        return
    end
    if isempty(hyphen) % single residue number
        try_numeric = str2double(range);
        if isnan(try_numeric) % this should be a residue type
            rtpoi = rtpoi + 1;
            residue_types{rtpoi} = range;
        else
            rpoi = rpoi + 1;
            residues(rpoi) = try_numeric;
        end
    else
        try_num_start = str2double(range(1:hyphen-1));
        try_num_end = str2double(range(hyphen+1:end));
        if isnan(try_num_start) || isnan(try_num_end)
            % should not happen, but allows for residue type names that do
            % contain a hyphen
            rtpoi = rtpoi + 1;
            residue_types{rtpoi} = range;
        else % residue number range
            rrange = str2double(range(1:hyphen-1)):str2double(range(hyphen+1:end));
            residues(rpoi+1:rpoi+length(rrange)) = rrange;
            rpoi = rpoi + length(rrange);
        end
    end
end
residues = residues(1:rpoi);
residue_types = residue_types(1:rtpoi);

% remove part of address that was already processed
if ~isempty(ei)
    address = address(ei+1:end);
end

% separate atoms part and make atom type cellstring
ei = strfind(address,':');
if ~isempty(ei)
    atomstr = address(1:ei-1);
end
delimited = strfind(atomstr,',');
atom_fields = cell(1,length(delimited)+1);
delimited = [0 delimited length(atomstr)+1];
for k = 1:length(atom_fields)
    atom_fields{k} = atomstr(delimited(k)+1:delimited(k+1)-1);
end

% remove part of address that was already processed
if ~isempty(ei)
    address = address(ei+1:end);
end

% make location string
locations = ' '; % the empty location tag is selected by default
if ~isempty(address)
    locations = ''; % if locations are specified,the empty location tag is not selected by default
    delimited = strfind(address,',');
    delimited = [0 delimited length(address)+1];
    for k = 1:length(delimited)-1
        locations = [locations address(delimited(k)+1:delimited(k+1)-1)]; %#ok<AGROW>
    end
end

exceptions = exceptions(1:warnings);

profile viewer