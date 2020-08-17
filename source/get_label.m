function [argsout,entity,exceptions] = get_label(entity,label,attributes,address,force)
%
% GET_LABEL Retrieves attributes of a label, computes the label if required
%
%   [argout,entity,exceptions] = GET_LABEL(entity,label,attributes)
%   Provides attribute values and possibly exceptions for sites selected
%   in an entity at residue level
%
%   [argout,entity,exceptions] = GET_LABEL(entity,label,attributes,address)
%   Provides attribute values and possibly exceptions for sites selected
%   by a residue address
%
%   [argout,emtity,exceptions] = GET_LABEL(entity,label,attributes,address,force)
%   if force = true, recomputation is performed even if the label had been
%   computed before, defaults to false
%
% INPUT
% entity       entity in an MMMx format, must be provided
% label        three-letter code or synonym for a label, either MMMx must
%              have a rotamer library for this label or it must be
%              atom.<atname>, where <atname> is a valid atom name in the
%              residue, such as 'CA' or 'Cu'
% attributes   residue attribute to be retrieved, string, defaults to
%              'info', either a string or a cell array of strings for
%              multiple simultaneous requests (saves time on selection)
%              attribute    output                           Type
%              ------------------------------------------------------------
%              info         rotamer library information      struct
%               .tlc        MMMx internal three-letter code  string
%                           for an atom site, ATM.<atname>
%               .class      label class, e.g., 'nitroxide'   string
%               .smiles     SMILES for label structure       string
%               .f-factor   forgive/attraction enhancement   (1,2) double
%               .ref_pop    reference populations for R      (R.1) double
%                           rotamers
%               .site       site address                     string
%              part_fun     partition function               double
%              torsion      T sidechain torsions for R       (R,T) double
%                           rotamers
%              orientations Euler angles of molecular frame  (R,3) double  
%              populations  populations for R rotamers       (R,1) double
%              positions    positions for R rotamers         (R,3) double
%              affine       affine transformation matrix     (4,4) double
%                           that transforms label to site
%              coor         full coordinates of R rotamers   (1,R) cell 
%              potentials   attachment energies (J/mol)      (R,1) double
%              numbers      numbers of the rotamers in the   (R,1) int
%                           library
% address      MMMx address, 'selected' refers to the current selection
%              defaults to 'selected'
% force        forces recomputation of the label rotamers if true, defaults
%              to false
%
% OUTPUT
% argsout      cell array of outputs, one cell for each site, see above
%              for a request with multiple attributes, cell array of cell
%              arrays, with outer level corresponding to attributes
% entity       input entity augmented by label information
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output
exceptions = {[]};
warnings = 0;

% reform single string input to cell aaray with on element
if ~iscell(attributes)
    attributes = {attributes};
end
num_attributes = length(attributes);
argsout = cell(1,num_attributes);

% convert label name into standard name

tlc = label_by_synonym(label);

% if there is no such label, check whether an tom is address
if isempty(tlc)
    if length(label) >= 5 && strcmpi(label(1:5),'atom.')
        atom_site = true;
        atom_tag = upper(label(6:end));
    else
        argsout = {};
        exceptions = {MException('get_label:unknown_label', 'Label %s is not known to MMMx',label)};
        return
    end
else
    atom_site = false;
    lib = load(sprintf('rotamer_library_%s.mat',tlc));
end


% set default arguments
if ~exist('address','var') || isempty(address)
    address = 'selected';
end
if ~exist('attributes','var') || isempty(attributes)
    attributes = {'info'};
end
if ~exist('force','var') || isempty(force)
    force = false;
end

% select the objects by the provided address
entity = select(entity,address,true);

index_vector = zeros(1,5,'uint16');
conformers = entity.selected; % selected conformers
% scan entity for selected locations
chains = fieldnames(entity);
for kc = 1:length(chains)
    chain = chains{kc};
    if isstrprop(chain(1),'upper') % chain fields start with a capital
        index_vector(1) =  entity.(chain).index;
        residues = fieldnames(entity.(chain));
        for kr = 1:length(residues) % expand over all residues
            residue = residues{kr};
            if strcmp(residue(1),'R') % these are residue fields
                index_vector(2) =  entity.(chain).(residue).index;
                if entity.(chain).(residue).selected % report back
                    % for an atom site, generate rotamer library equivalent
                    if atom_site
                        [lib,exceptions,warnings] = get_atom_rot_lib(entity,...
                            chain,residue,atom_tag,index_vector,conformers,exceptions,warnings);
                        if isempty(lib)
                            warnings = warnings + 1;
                            exceptions{warnings} = MException('get_label:unknown_atom_site',...
                                'Atom site %s does not exist',sprintf('(%s)%s.%s',chain,...
                                residue(2:end),atom_tag));
                            continue
                        end
                    end
                    for kattr = 1:num_attributes
                        attribute = attributes{kattr};
                        outputs = 0; % counter for the number of outputs
                        argout = cell(1,65000);
                        switch attribute
                            case 'info' % depends only on rotamer library
                                outputs = outputs + 1;
                                conformer_string = sprintf('%i',conformers(1));
                                for kconf = 2:length(conformers)
                                    conformer_string = sprintf('%s,%i',conformer_string,conformers(kconf));
                                end
                                info.tlc = lib.rot_lib.tlc;
                                info.class = lib.rot_lib.class;
                                info.smiles = lib.rot_lib.SMILES;
                                info.f_factor = lib.rot_lib.f_factor;
                                info.ref_pop = lib.rot_lib.populations;
                                info.site = sprintf('{%s}(%s)%s',conformer_string,chain,residue(2:end));
                                argout{outputs} = info;
                            case 'torsion' % depends only on rotamer library
                                outputs = outputs + 1;
                                if atom_site
                                    argout{outputs} = lib.torsion;
                                else
                                    torsion = lib.rot_lib.rotamers(1).torsion;
                                    torsion = zeros(1,length(torsion));
                                    for kt = 1:length(torsion)
                                        torsion(kt,:) = lib.rot_lib.rotamers(kt).torsion;
                                    end
                                    argout{outputs} = torsion;
                                end
                            case 'orientations' % depends on attachment
                                if atom_site
                                    outputs = outputs + 1;
                                    argout{outputs} = lib.orientations;
                                else
                                    for kconf = 1:length(conformers)
                                        entity = label_if_required(entity,chain,...
                                            residue,conformers(kconf),tlc,lib,force);
                                        outputs = outputs + 1;
                                        argout{outputs} = ...
                                            entity.(chain).(residue).labels.(tlc).orientations{conformers(kconf)};
                                    end
                                end
                            case 'populations' % depends on attachment
                                if atom_site
                                    outputs = outputs + 1;
                                    argout{outputs} = lib.populations;
                                else
                                    for kconf = 1:length(conformers)
                                        entity = label_if_required(entity,chain,...
                                            residue,conformers(kconf),tlc,lib,force);
                                        outputs = outputs + 1;
                                        argout{outputs} = ...
                                            entity.(chain).(residue).labels.(tlc).populations{conformers(kconf)};
                                    end
                                end
                            case 'positions' % depends on attachment
                                if atom_site
                                    outputs = outputs + 1;
                                    argout{outputs} = lib.positions;
                                else
                                    for kconf = 1:length(conformers)
                                        entity = label_if_required(entity,chain,...
                                            residue,conformers(kconf),tlc,lib,force);
                                        outputs = outputs + 1;
                                        argout{outputs} = ...
                                            entity.(chain).(residue).labels.(tlc).positions{conformers(kconf)};
                                    end
                                end
                            case 'affine' % depends on attachment
                                if atom_site
                                    outputs = outputs + 1;
                                    argout{outputs} = lib.affine;
                                else
                                    for kconf = 1:length(conformers)
                                        entity = label_if_required(entity,chain,...
                                            residue,conformers(kconf),tlc,lib,force);
                                        outputs = outputs + 1;
                                        argout{outputs} = ...
                                            entity.(chain).(residue).labels.(tlc).affine{conformers(kconf)};
                                    end
                                end
                            case 'part_fun' % depends on attachment
                                if atom_site
                                    outputs = outputs + 1;
                                    argout{outputs} = 1;
                                else
                                    for kconf = 1:length(conformers)
                                        entity = label_if_required(entity,chain,...
                                            residue,conformers(kconf),tlc,lib,force);
                                        outputs = outputs + 1;
                                        argout{outputs} = ...
                                            entity.(chain).(residue).labels.(tlc).part_fun{conformers(kconf)};
                                    end
                                end
                            case 'coor' % depends on attachment
                                if atom_site
                                    outputs = outputs + 1;
                                    argout{outputs} = lib.positions;
                                else
                                    for kconf = 1:length(conformers)
                                        entity = label_if_required(entity,chain,...
                                            residue,conformers(kconf),tlc,lib,force);
                                        outputs = outputs + 1;
                                        argout{outputs} = ...
                                            entity.(chain).(residue).labels.(tlc).coor{conformers(kconf)};
                                    end
                                end
                            case 'numbers' % depends on attachment
                                if atom_site
                                    outputs = outputs + 1;
                                    argout{outputs} = 1;
                                else
                                    for kconf = 1:length(conformers)
                                        entity = label_if_required(entity,chain,...
                                            residue,conformers(kconf),tlc,lib,force);
                                        outputs = outputs + 1;
                                        argout{outputs} = ...
                                            entity.(chain).(residue).labels.(tlc).numbers{conformers(kconf)};
                                    end
                                end
                            case 'potentials' % depends on attachment
                                if atom_site
                                    outputs = outputs + 1;
                                    argout{outputs} = 1;
                                else
                                    for kconf = 1:length(conformers)
                                        entity = label_if_required(entity,chain,...
                                            residue,conformers(kconf),tlc,lib,force);
                                        outputs = outputs + 1;
                                        argout{outputs} = ...
                                            entity.(chain).(residue).labels.(tlc).potentials{conformers(kconf)};
                                    end
                                end
                            otherwise
                                argsout = {};
                                exceptions = {MException('get_label:unsupported_attribute', 'Attribute %s not supported',attribute)};
                                return
                        end
                        argsout{kattr} = argout(1:outputs);
                    end
                end
            end
        end
    end
end

if num_attributes == 1
    argsout = argsout{1};
end

function [lib,exceptions,warnings] = get_atom_rot_lib(entity,chain,resname,...
    atom_tag,index_vector,conformers,exceptions,warnings)
% Makes a minimal pseudo-rotamer library for an atom site

residue = entity.(chain).(resname);

lib = '';

[m,~] = size(entity.index_array);
all_atom_indices = 1:m;

location_vector = 1:length(residue.locations);
populations0 = ones(length(location_vector),1);
% rotamers overrule locations
if length(residue.populations) > 1
    % all rotamers are selected
    location_vector = 1:length(residue.populations);
    populations0 = residue.populations;
end

positions = zeros(length(conformers)*length(location_vector),3);
populations = zeros(length(conformers)*length(location_vector),1);
outputs = 0;

% find all atom locations and the corresponding populations
if isfield(residue,atom_tag)
    index_vector(3) =  residue.(atom_tag).index;
    lib.rot_lib.SMILES = residue.(atom_tag).element;
    for kconf = 1:length(conformers)
        index_vector(4) = conformers(kconf);
        for kl = location_vector
            index_vector(5) = kl;
            [~,atom_indices] = ismember(entity.index_array,index_vector,'rows');
            atom_index = all_atom_indices(atom_indices~=0);
            if ~isempty(entity.xyz(atom_index,:))
                outputs = outputs + 1;
                positions(outputs,:) = entity.xyz(atom_index,:);
                populations(outputs) = populations0(kl)*entity.populations(kconf)*double(entity.occupancies(atom_index))/100;
            else
                warnings = warnings + 1;
                exceptions{warnings} = MException('get_label:unknown_atom_location',...
                    'Atom location %s does not exist',sprintf('{%i}(%s)%s.%s',conformers{kconf},chain,...
                    resname(2:end),atom_tag));
            end
        end
    end
end

% if no location of this atom was found, return empty pseudo-library
if outputs == 0
    lib = [];
    return
end

% set the pseduo-library entries
lib.populations = populations(1:outputs);
lib.rot_lib.populations = populations(1:outputs);
lib.positions = positions(1:outputs,:);
lib.orientations = [];
lib.torsion = [];
lib.rot_lib.class = 'atom';
lib.rot_lib.f_factor = [1,1];
lib.rot_lib.tlc = sprintf('ATM.%s',atom_tag);
lib.affine = eye(4);

function entity = label_if_required(entity,chain,residue,conformer,tlc,lib,force)
% Checks, whether rotamers were already computed for this label at this
% residue and computes them if nor or if forced

% check whether we need to compute or already have the result
compute = true;
if ~force && isfield(entity.(chain).(residue),'labels')
    if isfield(entity.(chain).(residue).labels,tlc)...
            && length(entity.(chain).(residue).labels.(tlc).affine) >= conformer...
            && ~isempty(entity.(chain).(residue).labels.(tlc).affine{conformer})
        compute = false;
    end
end

if compute
    site.conformer = conformer;
    site.chain = chain;
    site.residue = residue;
    rotamers = get_rotamers_for_site(entity,site,lib.rot_lib);
    entity.(chain).(residue).labels.(tlc).positions{conformer} = rotamers.positions;
    entity.(chain).(residue).labels.(tlc).populations{conformer} = rotamers.populations;
    entity.(chain).(residue).labels.(tlc).orientations{conformer} = rotamers.orientations;
    entity.(chain).(residue).labels.(tlc).affine{conformer} = rotamers.affine;
    entity.(chain).(residue).labels.(tlc).part_fun{conformer} = rotamers.part_fun;
    entity.(chain).(residue).labels.(tlc).coor{conformer} = rotamers.coor;
    entity.(chain).(residue).labels.(tlc).numbers{conformer} = rotamers.numbers;
    entity.(chain).(residue).labels.(tlc).potentials{conformer} = rotamers.potentials;
end