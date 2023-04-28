function [pairs,coulomb,sequence,resnums] = coulomb_pairs(entity,address,aa1,aa2,pH,I)
%
% COULOMB_PAIRS Returns ensemble_averaged Coulomb interaction for a pair of residue types
%               pH and ionic strength can be considered
%
%   [pairs,coulomb,sequence,resnums] = coulomb_pairs(entity,address,aa1,aa2)
%
%   [pairs,coulomb,sequence,resnums] = coulomb_pairs(entity,address,aa1,aa2,pH)
%
%   [pairs,coulomb,sequence,resnums] = coulomb_pairs(entity,address,aa1,aa2,pH,I)
%
%
% INPUT
% entity       entity in an MMMx format, must be provided
% address      MMMx address, defaults to chain 'A'
% aa1          type of first amino acid, defaults to 'Arg', case
%              insensitive
% aa2          type of second amino acid, defaults to 'Glu', case
%              insensitive
% pH           pH value, defaults to 7
% I            ionic strength, defaults to 150 mmol/L
%
% OUTPUT
% pairs        list of residue number pairs matching aa1 and aa2
% coulomb      ensemble averaged Coulomb interaction in temperature units [K],
%              negative for attraction, positive for repulsion
% sequence     amino acid sequence of the selected chain or residue range
% resnums      vector of residue numbers corresponding to sequence

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

epsilon0 = 8.8541878128e-12; % vacuum permittivity
epsilon = 80; % relative permittivity of water
kB = 1.380649e-23; % Boltzmann constant
elmc = 1.60217663e-19; % elementary charge

% default selection is chain A
if ~exist('address','var') || isempty(address)
    address = '(A)';
end

% default amino acid type 1 is 'Arg'
if ~exist('aa1','var') || isempty(aa1)
    aa1 = 'Arg';
end

% default amino acid type 2 is 'Glu'
if ~exist('aa2','var') || isempty(aa2)
    aa2 = 'Glu';
end

% default pH is 7
if ~exist('pH','var') || isempty(pH)
    pH = 7;
end

if ~exist('I','var') || isempty(I)
    I = 0.150;
end
lambda_D = debye_length(I);

pop = entity.populations;
C = length(entity.populations); % number of conformers

[selected_chain,range] = split_chain_range(address);
sequence = '';
resnums = [];

chains = fieldnames(entity);
for kc = 1:length(chains)
    chain = chains{kc};
    if isstrprop(chain(1),'upper') && strcmpi(chain(1),selected_chain) % chain fields start with a capital
        residues = fieldnames(entity.(chain));
        resnums = zeros(length(residues),1);
        coulomb = zeros(length(residues),1);
        pairs = zeros(length(residues),2);
        respoi = 0;
        pairpoi = 0;
        for kr = 1:length(residues) % expand over all residues
            residue = residues{kr};
            resnum = str2double(residue(2:end));
            if strcmp(residue(1),'R') && resnum >= range(1) && resnum <= range(2) % these are residue fields
                resname = entity.(chain).(residue).name;
                slc = tlc2slc(resname);
                sequence = strcat(sequence,slc);
                respoi = respoi + 1;
                resnums(respoi) = resnum; 
                if strcmpi(resname,aa1)
                    switch resname
                        case 'ASP'
                            q1 = get_charge(3.90,pH) - 1;
                            indices1 = zeros(2,C);
                            indices1(1,:) = entity.(chain).(residue).OD1.tab_indices;
                            indices1(2,:) = entity.(chain).(residue).OD2.tab_indices;
                        case 'GLU'
                            q1 = get_charge(4.07,pH) - 1;
                            indices1 = zeros(2,C);
                            indices1(1,:) = entity.(chain).(residue).OE1.tab_indices;
                            indices1(2,:) = entity.(chain).(residue).OE2.tab_indices;
                        case 'HIS'
                            q1 = get_charge(6.04,pH);
                            indices1 = zeros(2,C);
                            indices1(1,:) = entity.(chain).(residue).ND1.tab_indices;
                            indices1(2,:) = entity.(chain).(residue).NE2.tab_indices;
                        case 'LYS'
                            q1 = get_charge(10.54,pH);
                            indices1 = entity.(chain).(residue).NZ.tab_indices;
                        case 'ARG'
                            q1 = get_charge(12.48,pH);
                            indices1 = entity.(chain).(residue).CZ.tab_indices;
                        otherwise
                            q1 = 0;
                            indices1 = [];
                    end
                    for kr2 = 1:length(residues) % expand over all residues
                        residue2 = residues{kr2};
                        resnum2 = str2double(residue2(2:end));
                        if strcmp(residue2(1),'R') && resnum2 >= range(1) && resnum2 <= range(2) && kr2 ~= kr % these are residue fields
                            resname2 = entity.(chain).(residue2).name;
                            if strcmpi(resname2,aa2)
                                pairpoi = pairpoi + 1;
                                pairs(pairpoi,1) = resnum;
                                pairs(pairpoi,2) = resnum2;
                                switch resname2
                                    case 'ASP'
                                        q2 = get_charge(3.90,pH) - 1;
                                        indices2 = zeros(2,C);
                                        indices2(1,:) = entity.(chain).(residue2).OD1.tab_indices;
                                        indices2(2,:) = entity.(chain).(residue2).OD2.tab_indices;
                                    case 'GLU'
                                        q2 = get_charge(4.07,pH) - 1;
                                        indices2 = zeros(2,C);
                                        indices2(1,:) = entity.(chain).(residue2).OE1.tab_indices;
                                        indices2(2,:) = entity.(chain).(residue2).OE2.tab_indices;
                                    case 'HIS'
                                        q2 = get_charge(6.04,pH);
                                        indices2 = zeros(2,C);
                                        indices2(1,:) = entity.(chain).(residue2).ND1.tab_indices;
                                        indices2(2,:) = entity.(chain).(residue2).NE2.tab_indices;
                                    case 'LYS'
                                        q2 = get_charge(10.54,pH);
                                        indices2 = entity.(chain).(residue2).NZ.tab_indices;
                                    case 'ARG'
                                        q2 = get_charge(12.48,pH);
                                        indices2 = entity.(chain).(residue2).CZ.tab_indices;
                                    otherwise
                                        q2 = 0;
                                        indices2 = [];
                                end
                                if q1*q2 ~= 0
                                    ECoulomb = 0;
                                    for c = 1:C % loop over conformers
                                        cindices1 = indices1(:,c); % atom indices for this residue in this conformer
                                        coor1 = entity.xyz(cindices1,:); % atom coordinates for this residue in this conformer
                                        coor1 = mean(coor1,1); % mean coordinate
                                        cindices2 = indices2(:,c); % atom indices for this residue in this conformer
                                        coor2 = entity.xyz(cindices2,:); % atom coordinates for this residue in this conformer
                                        coor2 = mean(coor2,1); % mean coordinate
                                        r = norm(coor1-coor2);
                                        ECoulomb = ECoulomb + pop(c)*elmc^2*q1*q2*exp(-r/lambda_D)/(4*pi*epsilon0*epsilon*1e-10*r);
                                    end
                                    coulomb(pairpoi) = ECoulomb/kB;
                                end
                            end
                        end
                    end                    
                end
            end % is a residue            
        end % residue fields loop
    end % is a chain
end
coulomb = coulomb(1:pairpoi);
pairs = pairs(1:pairpoi,:);
resnums = resnums(1:respoi);

function charge = get_charge(pKa,pH)

charge = exp(-log(10)*(pH-pKa))/(1+exp(-log(10)*(pH-pKa)));

function [chain,range] = split_chain_range(address)

range = [1,1e6];

poia = strfind(address,'(');
poie = strfind(address,')');
if isempty(poie)
    chain = '';
else
    chain = address(poia+1:poie-1);
    address = address(poie+1:end);
end
if strcmp(chain,'*')
    return
end
residues = split(address,'-');
if ~isempty(residues{1})
    range(1) = str2double(residues{1});
end
if length(residues)>1 && ~isempty(residues{2})
    range(2) = str2double(residues{2});
end
if length(range) == 1
    range(2) = range(1);
end
