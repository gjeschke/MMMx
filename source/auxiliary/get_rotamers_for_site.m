function rotamers = get_rotamers_for_site(entity,site,rot_lib)
%
% GET_ROTAMERS_FOR_SITE computes label rotamers at a single site
%        
%
%   [rotamers,exceptions] = GET_ROTAMERS_FOR_SITE(entity,site,rot_lib)
%   Rotamer information for a site in entity for a label specified by a
%   rotamer library rot_lib
%
% INPUT
% entity     entity in MMMx:atomic representation, must contain all context
% site       labeling site specification, struct with fields
%            .conformer conformer number, defaults to 1
%            .chain     chain identifier
%            .residue   residue number or identifier
% rot_lib    rotamer library for the label in MMMx format
%
% OUTPUT
% rotamers   structure with information on the attached rotamers, rotamers
%            are sorted by descending population
%            .positions      (R,3) double coordinates of R rotamers making
%                            up 99.5% of total rotamer population
%            .populations    (R,1) double populations of R rotamers
%            .orientations   (R,3) Euler angles relating the label frame to
%                            the entity frame for R rotamers
%            .potentials     (R,1) double entity interaction energies
%                            (J/mol) for R rotamers
%            .numbers        (R,1) int numbers of the R rotamers in the
%                            library
%            .coor           (1,R) cell with full coordinates of R rotamers
%            .part_fun       double partition function for attachment
%            .affine         affine transformation from standard frame to
%                            local site frame 
% exceptions   cell vector of MException objects if something went wrong, 
%              defaults to one cell holding an empty array
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2010-2020: Yevhen O. Polyhach, Gunnar Jeschke

% constants
gas_un=8.314472;    % universal gas constant in CI (J/(mol*K)       
T = 298; % ambient temperature [K]

threshold=0.005; % 0.5% of population is neglected

% Initialize empty output
rotamers = [];

% set defaults for missing input
if ~isfield(site,'conformer') || isempty(site.conformer)
    site.conformer = 1;
end

% reformat residue number to residue tag, if necessary
if ~ischar(site.residue)
    site.residue = sprintf('R%i',site.residue);
end

% run options for interaction energy computation
options.ff = rot_lib.ff;
options.f_factor = rot_lib.f_factor;

% make context coordinates including atomic numbers
% select all atoms of all redidues residues in the whole conformer
entity = select(entity,sprintf('{%i}(*)*.*',site.conformer));
% unselect the site
entity = select(entity,sprintf('{%i}(%s)%s.*',site.conformer,site.chain,site.residue(2:end)),false,true);
% argout = get_atom(entity,'xyz');
% arrayout = argout(~cellfun('isempty',argout));
% outputs = length(arrayout);
% xyz = cell2mat(arrayout);
% xyz = reshape(xyz,3,outputs).';
% argout = get_atom(entity,'element');
% elements = double(cell2mat(argout).');
% context = [elements xyz];

argouts = get_atom(entity,'ecoor');
context = argouts{1}(:,1:4);

% get standard frame coordinates
address = sprintf('{%i}(%s)%s.%s',site.conformer,site.chain,site.residue(2:end),...
    rot_lib.std_frame_atoms{1});
orig = get_atom(entity,'xyz',address);
orig = orig{1};
address = sprintf('{%i}(%s)%s.%s',site.conformer,site.chain,site.residue(2:end),...
    rot_lib.std_frame_atoms{2});
xax = get_atom(entity,'xyz',address);
xax = xax{1};
address = sprintf('{%i}(%s)%s.%s',site.conformer,site.chain,site.residue(2:end),...
    rot_lib.std_frame_atoms{3});
ypax = get_atom(entity,'xyz',address);
ypax = ypax{1};
x = xax - orig; 
x = x/norm(x);    
yp = ypax - orig; 
yp = yp/norm(yp);
z = cross_rowvec(x,yp); 
z = z/norm(z);
y = cross_rowvec(z,x);  
Rp = [x;y;z]; % rotation matrix for conversion to standard frame
% generate an affine transformation matrix and store it
transmat1 = eye(4);
transmat1(1:3,4) = orig.';
transmat2 = eye(4);
transmat2(1:3,1:3) = Rp.';
transmat = transmat1*transmat2;
rotamers.affine = transmat;

% number of rotamers
R = length(rot_lib.rotamers);
rotamers.orientations = zeros(R,3);
rotamers.positions = zeros(R,3);
rotamers.potentials = zeros(R,1);
rotamers.coor = cell(1,R);
% number of atoms
[A,~] = size(rot_lib.rotamers(1).coor);
% first sidechain atom
first_side_chain_atom = rot_lib.side_chain;
% atomic number vector for sidechain
sc_elements = rot_lib.elements(first_side_chain_atom:end);
% transform rotamers and determine potentials
for r = 1:R
    coor = [rot_lib.rotamers(r).coor.';ones(1,A)];
    coor = transmat*coor;
    coor = coor(1:3,:).';
    rotamers.coor{r} = coor;
    molframe = coor(rot_lib.mol_frame,:);
    x = molframe(2,:) - molframe(1,:);
    x = x/norm(x);
    yp = molframe(3,:) - molframe(1,:);
    yp = yp/norm(yp);
    z = cross_rowvec(x,yp);
    z = z/norm(z);
    y = cross_rowvec(z,x);
    rotamers.orientations(r,:) = DCM2Euler([x;y;z]);
    side_chain = [sc_elements coor(first_side_chain_atom:end,:)];
    rotamers.potentials(r) = energy_LJ(side_chain,context,options);
end
populations = exp(-rotamers.potentials/(gas_un*T));
rotamers.populations = populations.*rot_lib.populations;
% partition function for attachment
rotamers.part_fun = sum(rotamers.populations)/sum(rot_lib.populations);
% renormalize populations
rotamers.populations = rotamers.populations/sum(rotamers.populations);
% sort rotamers by descending population
[rotamers.populations,sorting] = sort(rotamers.populations,'descend');
rotamers.coor = rotamers.coor(sorting);
% remove the least populated rotamers above a coverage threshold
decider = cumsum(rotamers.populations);
[~,significant_rotamers] = min(abs(decider-1+threshold));
if decider(significant_rotamers) < 1- threshold
    significant_rotamers = significant_rotamers + 1;
end
rotamers.populations = rotamers.populations(1:significant_rotamers);
rotamers.populations = rotamers.populations/sum(rotamers.populations);
rotamers.coor = rotamers.coor(1:significant_rotamers);
rotamers.numbers = sorting(1:significant_rotamers);

% compute label positions
position_indices =rot_lib.position(:,1);
position_weights =rot_lib.position(:,2);
% renormalize position weights
position_weights = position_weights/sum(position_weights);
% positions for all rotamers
for r = 1:significant_rotamers
    coor = rotamers.coor{r}; % r-th rotamer
    position = coor(position_indices,:);
    position = position_weights'*position; % mean position coordinates for the r-th rotamer
    rotamers.positions(r,:) = position; 
end

rotamers.positions = rotamers.positions(1:significant_rotamers,:);

