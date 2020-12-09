function [entity,exceptions] = select(entity,address,overwrite,unselect)
%
% SELECT Select objects in an entity
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


exceptions{100} = [];
warnings = 0; % counter for warnings

% initialize overwrite, if it does not exist
if ~exist('overwrite','var') || isempty(overwrite)
    overwrite = false;
end

% initialize invert, if it does not exist
if ~exist('unselect','var') || isempty(unselect)
    unselect = false;
elseif unselect
    overwrite = false;
end

% treat address `selected`
if strcmp(address,'selected') 
    if ~unselect
        exceptions{1} = [];
        return
    else
        overwrite = true;
    end
end

if overwrite % delete existing selection
    entity.selected = 1;
    entity.water_selected = false;
    chains = fieldnames(entity);
    for kc = 1:length(chains)
        chain = chains{kc};
        if isstrprop(chain(1),'upper') % chain fields start with a capital
            entity.(chain).selected = 0;
            residues = fieldnames(entity.(chain));
            for kr = 1:length(residues)
                residue = residues{kr};
                if strcmp(residue(1),'R') % these are residue fields
                    entity.(chain).(residue).selected = 0;
                    entity.(chain).(residue).selected_rotamers = 1;
                    atoms = fieldnames(entity.(chain).(residue));
                    for ka = 1:length(atoms)
                        atom = atoms{ka};
                        if isstrprop(atom(1),'upper') % these are atom fields
                            entity.(chain).(residue).(atom).selected = 0;
                            entity.(chain).(residue).(atom).selected_locations = 1;
                        end
                    end
                end
            end
        end
    end
end

% we have now unselected the selection, if this was requested, and can return
if strcmp(address,'selected') 
    exceptions{1} = [];
    return
end

% treat water selection
if strcmpi(address,'water')
    entity.water_selected = true;
    exceptions = {};
    return
end

% extract conformer selection if any
si = strfind(address,'{');
ei = strfind(address,'}');
if ~isempty(si) && ~isempty(ei)
    confstr = address(si+1:ei-1);
    address = [address(1:si-1) address(ei+1:end)];
else
    confstr = '';
end
if strcmp(strtrim(confstr),'*')
    all_conformers = true;
    if unselect
        entity.selected = 1;
    else
        entity.selected =  1:length(entity.populations);
    end
else
    all_conformers = false;
end
% make conformer vector
if ~isempty(confstr) && ~all_conformers
    conformers = zeros(1,66000);
    delimited = strfind(confstr,',');
    delimited = [0 delimited length(confstr)+1];
    cpoi = 0;
    for k = 1:length(delimited)-1
        range = confstr(delimited(k)+1:delimited(k+1)-1);
        hyphen = strfind(range,'-');
        if length(hyphen) > 1
            exceptions = {MException('select:wrong_conformer_range', 'Conformer range contains more than one hyphen')};
            return
        end
        if isempty(hyphen) % single conformer
            cpoi = cpoi + 1;
            % process start/end syntax
            if strcmp(strtrim(range),'start')
                conformers(cpoi) = 1;
            elseif strcmp(strtrim(range),'end')
                conformers(cpoi) = length(entity.populations);
            else
                conformers(cpoi) = str2double(range);
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
            conformers(cpoi+1:cpoi+length(crange)) = crange;
            cpoi = cpoi + length(crange);
        end
    end
    if unselect
        old_conformers = entity.selected;
        new_conformers = old_conformers;
        ncpoi = 0;
        for kconf = 1:length(old_conformers)
            % if it does not match, it remains selected
            if min(abs(conformers(1:cpoi)-old_conformers(kconf))) > 0
                ncpoi = ncpoi + 1;
                new_conformers(ncpoi) = old_conformers(kconf);
            end
        end
        if ncpoi == 0
            new_conformers = 1; % minimal conformer selection
        else
            new_conformers = new_conformers(1:ncpoi); % the ones still selected
        end
        entity.selected = new_conformers;
    elseif overwrite
        entity.selected = conformers(1:cpoi);
    else
        is_added = true(1,length(entity.selected));
        for ocpoi = 1:length(entity.selected)
            for ncpoi = 1:cpoi
                if entity.selected(ocpoi) == conformers(ncpoi)
                    is_added(ocpoi) = false;
                end
            end
        end
        entity.selected = [entity.selected(is_added) conformers(1:cpoi)];
    end
end

level = 0; % initializes selection level to entity
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
    level = 1;
    address = address(ei+1:end);
end

% separate residues part and make residue vector and residue type
% cellstring
residue_numbers = zeros(1,66000);
rpoi = 0;
residue_types = cell(1,100);
rtpoi = 0;
irot = strfind(address,'|'); % check whether there is a rotamer specification
ei = strfind(address,'.');
if ~isempty(address) && isempty(ei)
    ei = length(address)+1;
end
rotamer_str = '';
erot = [];
if ~isempty(irot)
    rotamer_str = address(irot+1:ei-1);
    erot = ei;
    ei = irot;
end
if ~isempty(ei)
    resstr = address(1:ei-1);
    delimited = strfind(resstr,',');
    delimited = [0 delimited length(resstr)+1];
    for k = 1:length(delimited)-1
        range = resstr(delimited(k)+1:delimited(k+1)-1);
        hyphen = strfind(range,'-');
        if length(hyphen) > 1
            exceptions = {MException('select:wrong_residue_range', 'Residue range contains more than one hyphen')};
            return
        end
        if isempty(hyphen) % single residue number
            try_numeric = str2double(range);
            if isnan(try_numeric) % this should be a residue type
                rtpoi = rtpoi + 1;
                residue_types{rtpoi} = range;
            else
                rpoi = rpoi + 1;
                residue_numbers(rpoi) = try_numeric;
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
                residue_numbers(rpoi+1:rpoi+length(rrange)) = rrange;
                rpoi = rpoi + length(rrange);
            end
        end
    end
    residue_numbers = residue_numbers(1:rpoi);
    residue_types = residue_types(1:rtpoi);
    level = 2; % selection level residues
else
    residue_numbers = [];
    residue_types = {};
end

% remove part of address that was already processed
if ~isempty(ei)
    address = address(ei+1:end);
    erot = erot-ei;
end

% process rotamer selection if any
if strcmp(strtrim(rotamer_str),'*') % wildcard selects all rotamers
    all_rotamers = true;
else
    all_rotamers = false;
end
if ~isempty(rotamer_str) && ~all_rotamers
    rotamers = zeros(1,70000);
    rotpoi = 0;
    delimited = strfind(rotamer_str,',');
    delimited = [0 delimited length(rotamer_str)+1];
    for k = 1:length(delimited)-1
        range = rotamer_str(delimited(k)+1:delimited(k+1)-1);
        hyphen = strfind(range,'-');
        if length(hyphen) > 1
            exceptions = {MException('select:wrong_rotamer_range', 'Rotamer range contains more than one hyphen')};
            return
        end
        if isempty(hyphen) % single residue number
            try_numeric = str2double(range);
            if isnan(try_numeric) % wrong rotamer selection
                exceptions = {MException('select:wrong_rotamer_identifier', 'Rotamer identifier %s is not numeric', range)};
                return
            else
                rotpoi = rotpoi + 1;
                rotamers(rotpoi) = try_numeric;
            end
        else
            try_num_start = str2double(range(1:hyphen-1));
            try_num_end = str2double(range(hyphen+1:end));
            if isnan(try_num_start) || isnan(try_num_end)
                exceptions = {MException('select:wrong_rotamer_identifier', 'Rotamer identifier %s is not numeric', range)};
                return
            else % rotamer number range
                rotrange = try_num_start:try_num_end;
                rotamers(rotpoi+1:rotpoi+length(rotrange)) = rotrange;
                rotpoi = rotpoi + length(rotrange);
            end
        end
    end
    rotamers = rotamers(1:rotpoi);
else
    rotamers = 1; % if no rotamers are specified, rotamer 1 is selected by default
end
if ~isempty(erot)
    address = address(erot+1:end);
end
% separate atoms part and make atom type cellstring
ei = strfind(address,':');
if isempty(ei) && ~isempty(address)
    ei = length(address) + 1;
end
if ~isempty(ei)
    atomstr = address(1:ei-1);
    delimited = strfind(atomstr,',');
    atom_fields = cell(1,length(delimited)+4);
    delimited = [0 delimited length(atomstr)+1];
    apoi = 0;
    for k = 1:length(delimited)-1
        atom_type = atomstr(delimited(k)+1:delimited(k+1)-1);
        if strcmpi(atom_type,'backbone')
            atom_fields{apoi+1} = 'N';
            atom_fields{apoi+2} = 'CA';
            atom_fields{apoi+3} = 'C';
            atom_fields{apoi+4} = 'O';
            apoi = apoi + 4;
        else
            apoi = apoi + 1;
            atom_fields{apoi} = atom_type;
        end
    end
    atom_fields = atom_fields(1:apoi);
    % remove part of address that was already processed
    if ~isempty(ei)
        level = 3; % selection level atoms
        address = address(ei+1:end);
    end
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
    level = 4;
end

exceptions = exceptions(1:warnings);

if level == 0 % nothing was selected
    return
end

% selection is set in entity

chains = fieldnames(entity);
for kc = 1:length(chains)
    chain = chains{kc};
    if isstrprop(chain(1),'upper') % chain fields start with a capital
        if level == 1
            selected = entity.(chain).selected; % current selection state
            assignment = ~unselect; % select or unselect
            for ksel = 1:length(chain_fields)
                if strcmp(chain_fields{ksel},chain) ... % new selection matches?
                        || strcmp(chain_fields{ksel},'*') % wildcard selects all chains
                    selected = assignment;
                end
            end
            entity.(chain).selected = selected;
        else
            active = isempty(chain_fields{1});
            for ksel = 1:length(chain_fields)
                if strcmp(chain_fields{ksel},chain) ... % new selection matches?
                        || strcmp(chain_fields{ksel},'*') % wildcard selects all chains
                    active = true;
                end
            end
            if ~active
                continue;
            end
            residues = fieldnames(entity.(chain));
            for kr = 1:length(residues)
                residue = residues{kr};
                if strcmp(residue(1),'R') % these are residue fields
                    if all_rotamers % assign rotamer vector for wildcard
                        rotamers = 1:length(entity.(chain).(residue).populations);
                    end
                    location_set = entity.(chain).(residue).locations;
                    if level == 2
                        selected = entity.(chain).(residue).selected; % current selection state
                        was_selected = selected; % remember it for rotamer reassignment
                        assignment = ~unselect; % select or unselect
                        for ksel = 1:length(residue_numbers)
                            if residue_numbers(ksel) == str2double(residue(2:end)) % residue number match?
                                selected = assignment;
                            end
                        end
                        for ksel = 1:length(residue_types)
                            if strcmpi(residue_types(ksel),entity.(chain).(residue).name) ... % residue type match?
                                    || strcmp(residue_types{ksel},'*') % wildcard selects all residues
                                selected = assignment;
                            end
                        end
                        entity.(chain).(residue).selected = selected;
                        if selected && ~was_selected % rotamer assignment for new selections
                            entity.(chain).(residue).selected_rotamers...
                                = rotamers;
                        elseif ~selected % remove rotamer assignment if residue is not selected
                            entity.(chain).(residue).selected_rotamers...
                                = 1;
                        end
                    else
                        active = isempty(residue_numbers);
                        if active
                            active = isempty(residue_types{1});
                        end
                        for ksel = 1:length(residue_numbers)
                            if residue_numbers(ksel) == str2double(residue(2:end)) % residue number match?
                                active = true;
                            end
                        end
                        for ksel = 1:length(residue_types)
                            if strcmpi(residue_types(ksel),entity.(chain).(residue).name) ... % new selection matches?
                                    || strcmp(residue_types{ksel},'*') % wildcard selects all residues
                                active = true;
                            end
                        end
                        if ~active
                            continue;
                        end
                        % process location wildcard
                        if strcmp(strtrim(locations),'*')
                            res_locations = entity.(chain).(residue).locations;
                        else
                            res_locations = locations;
                        end
                        % translate location string to location index vector
                        location_vector = zeros(1,length(res_locations));
                        locpoi = 0;
                        for kloc = 1:length(res_locations)
                            loc = strfind(location_set,res_locations(kloc));
                            if ~isempty(loc)
                                locpoi = locpoi + 1;
                                location_vector(locpoi) = loc;
                            end
                        end
                        location_vector = location_vector(1:locpoi);
                        % rotamer selection overrules location selection,
                        % if both are present
                        if length(rotamers) > 1 || rotamers(1) ~= 1 
                            location_vector = rotamers;
                        end
                        % at least one location must be selected
                        no_location = false;
                        if isempty(location_vector)
                            location_vector = 1;
                            if level == 4
                                no_location = true;
                            end
                        end
                        atoms = fieldnames(entity.(chain).(residue));
                        for ka = 1:length(atoms)
                            atom = atoms{ka};
                            if isstrprop(atom(1),'upper') % these are atom fields
                                selected = entity.(chain).(residue).(atom).selected; % old selection state
                                was_selected = selected; % memorize for location assignment
                                assignment = ~unselect; % select or unselect?
                                for ksel = 1:length(atom_fields)
                                    if strcmp(atom_fields{ksel},atom) ... % atom type match?
                                            || strcmp(atom_fields{ksel},'*') % wildcard selects all atoms
                                        selected = assignment;
                                    end
                                end
                                entity.(chain).(residue).(atom).selected...
                                    = selected;
                                if no_location
                                    selected = false;
                                    entity.(chain).(residue).(atom).selected = false;
                                end
                                if selected && ~was_selected % loction assignment for new selection
                                    entity.(chain).(residue).(atom).selected_locations...
                                        = location_vector;
                                elseif ~selected % default location if atom is now unselected
                                    entity.(chain).(residue).(atom).selected_locations...
                                        = 1;
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
