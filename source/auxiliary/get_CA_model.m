function CA_model = get_CA_model(entity,conformer)
% CA_model = GET_CA_MODEL(entity,conformer)
%
% extracts Calpha coordinates, residue masses, assigned residue addresses, 
% Bfactors, and residue types of a Calpha trace model from an MMMx entity
%
% if argument modnum is empty, conformer 1 is addressed
%
% an empty output is returned if the entity has no peptide chains,
% if a Ca atom has alternate locations, an average coordinate is returned
% any residues other than amino acids are ignored
% if the structure does not contain hydrogen atoms, only heavy-atom masses
% are considered
%
% indices of the residues and Bfactors are returned for later reference
%
% Input:
%
% entity    entity in MMMx:atomic format
% conformer number of the conformer, for which the model is extracted,
%           defaults to 1
%
% Output:
%
% CA_model  struct with information on the Calpha trace model, fields
%           .coor       array(N,3) of Cartesian coordinates of Calpha atoms
%           .masses     array(N,1) of residue masses
%           .bfactors   array(N,1) with isotropic Bfactors, NaN if missing
%           .chains     array(N,1) of chain indices
%           .addresses  cell(N,1) of corresponding residue addresses

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% G. Jeschke, 2010-2021

% load the periodic system of elements and residue attributes, as we need 
% them for mass determination
chem = load('element_attributes.mat');
bio = load('monomer_attributes.mat');

if ~exist('conformer','var') || isempty(conformer)
    conformer = 1;
end

% initialize empty output
CA_model.coor = [];
CA_model.masses = [];
CA_model.bfactors = [];
CA_model.chains = [];
CA_model.addresses = {};

% initialize temporary output
max_nodes = 20000;
nodes = 0;
coor = zeros(max_nodes,3);
masses = zeros(max_nodes,1);
chain_indices = zeros(max_nodes,1);
bfactors = NaN(max_nodes,1);
addresses = cell(max_nodes,1);

% extract data from entity
chains = fieldnames(entity);
for kc = 1:length(chains) % loop over all chains
    chain = chains{kc};
        if isstrprop(chain(1),'upper') % chain fields start with a capital
            entity.(chain).selected = 0;
            residues = fieldnames(entity.(chain));
            for kr = 1:length(residues)
                residue = residues{kr};
                if strcmp(residue(1),'R') % these are residue fields
                    % consider a residue as peptide residue, if it has at
                    % least N, CA, and C atoms
                    if isfield(entity.(chain).(residue),'N') &&...
                            isfield(entity.(chain).(residue),'CA') &&...
                            isfield(entity.(chain).(residue),'C')
                        nodes = nodes + 1; % this is a new network node
                        chain_indices(nodes) = entity.(chain).index;
                        addresses{nodes} = sprintf('(%s)%s',chain,residue(2:end));
                        % select the table entry of the CA atom
                        tab_index = entity.(chain).(residue).CA.tab_indices(conformer);
                        % extract coordinate
                        coor(nodes,:) = entity.xyz(tab_index,:);
                        % extract B factor, if provided, else keep it NaN
                        if entity.(chain).(residue).CA.bfactor ~= 0
                            bfactors(nodes) = entity.(chain).(residue).CA.bfactor;
                        end
                        % determine residue mass
                        % if it is a known residue, use known residue mass
                        restype = tag2id(entity.(chain).(residue).name,upper(bio.monomers.aa_tags));
                        if ~isempty(restype)
                            mass = bio.monomers.aa.mass(restype);
                        else % otherwise determine mass from resolved atoms
                            mass = 0;
                            atoms = fieldnames(entity.(chain).(residue));
                            for ka = 1:length(atoms)
                                atom = atoms{ka};
                                if isstrprop(atom(1),'upper') % these are atom fields
                                    tab_index = entity.(chain).(residue).(atom).tab_indices(conformer);
                                    element = entity.elements(tab_index);
                                    mass = mass + chem.pse.mass(element);
                                end
                            end
                        end
                        masses(nodes) = mass;
                    end
                end
            end
        end
end

if nodes > 0
    % restrict arrays to actual number of amino acid residues
    CA_model.coor = coor(1:nodes,:);
    CA_model.masses = masses(1:nodes);
    CA_model.chains = chain_indices(1:nodes);
    CA_model.bfactors = bfactors(1:nodes);
    CA_model.addresses = addresses(1:nodes);
end

