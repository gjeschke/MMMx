function [pre_list,exceptions,entity] = pre_of_entity(entity,site_list,label,label_adr,td,taur,taui,R2_dia,larmor)
%
% PRE_OF_ENTITY Computation of paramagnetic relaxation enhancement for
%               selected NH protons of an entity
%
%   pre_list = PRE_OF_entity(entity,site_list,label_coor,taur,taui,R2_dia)
%   Given an entity, a list of sites, and a label type, the protein
%   rotational correlation time and label internal correlation time and the
%   nuclear transverse relaxation rate for the diamagnetically labeled
%   protein, intensity rations I_para/I_dia are estimated
%
%   see: Tesei et al. https://doi.org/10.1101/2020.08.09.243030 
%        for description of the approach
%
% INPUT
% entity    MMMx entity structure of a protein (ensemble)
% site_list list of the residues for which the PREs should be computed,
%           vector (1,n_sites) of struct
%           .chain    char, chain identifier
%           .residue  int, residue number
% label     label type, defaults to 'mtsl'
% label_adr label address without conformer specification. e.g. (A)316
% td        double, total INEPT time in HSQC [s]
% taur      double, rotational correlation time of the protein [s]
% taui      double, correlation time of internal label motion [s]
% R2_dia    nuclear transverse relataion rate(s), can be either one double
%           value applying to all nuclei or a vector (1,N) of values for
%           individual nuclei [s^{-1)]
% larmor    larmor frequency [MHz]
%
% OUTPUT
% pre_list      struct (n_sites,1) of
%               .chain      chain identifier
%               .residue    residue number
%               .Gamma2     PRE rate
%               .pre        PRE ratio
% exceptions    cell vector of MException objects if something went wrong, 
%               defaults to one cell holding an empty array
% entity        entity, which may be augmented by a label
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2010-2020: Yevhen Polyhach, Gunnar Jeschke

NHbond = 0.98; % length of N-H bond in Angstroem for proton construction

if isempty(label)
    label = 'mtsl';
end

% initialize output
n_sites = length(site_list);
pre_list(n_sites).chain = '';
pre_list(n_sites).residue = 0;
pre_list(n_sites).Gamma2 = [];
pre_list(n_sites).pre = [];
exceptions = {[]};
warnings = 0;

% pres = zeros(n_sites,1);
Gamma2 = zeros(n_sites,1);

for c = 1:length(entity.populations) % loop over all conformers
    % generate coordinate array from site list
    coor = zeros(n_sites,3);
    n_H = 0; % counter for protons that could be found, for instance, proline does not have one
    for site = 1:n_sites
        if ~isfield(entity,site_list(site).chain) % if a wrong chain is addressed, skip entry and raise warning
            warnings = warnings + 1;
            exceptions{warnings} = MException('pre_of_entity:unknown_chain',...
                'Chain %s does not exist',site_list(site).chain);
            continue
        end
        resfield = sprintf('R%i',site_list(site).residue);
        if ~isfield(entity.(site_list(site).chain),resfield) % if a wrong residue is addressed, skip entry and raise warning
            warnings = warnings + 1;
            exceptions{warnings} = MException('pre_of_entity:unknown_residue',...
                'Residue %i does not exist in chain %s',site_list(site).residue,site_list(site).chain);
            continue
        end
        % check, if the structure contains the coordinate of the NH proton
        if isfield(entity.(site_list(site).chain).(resfield),'H')
            H_index = entity.(site_list(site).chain).(resfield).H.tab_indices(c);
            n_H = n_H + 1;
            coor(n_H,:) = entity.xyz(H_index,:);
            pre_list(n_H).chain = site_list(site).chain;
            pre_list(n_H).residue = site_list(site).residue;
        else % otherwise, try to construct it
            resfield_prev = sprintf('R%i',site_list(site).residue - 1);
            % if proton coordinate does not exist and previous residue in the
            % chain does not exist, skip and raise warning
            if ~isfield(entity.(site_list(site).chain),resfield_prev)
                warnings = warnings + 1;
                exceptions{warnings} = MException('pre_of_entity:no_proton',...
                    'For residue %i in chain %s proton is missing and cannot be constructed',site_list(site).residue,site_list(site).chain);
                continue
            else
                failed = false;
                if isfield(entity.(site_list(site).chain).(resfield),'N')
                    N_index = entity.(site_list(site).chain).(resfield).N.tab_indices(c);
                else
                    failed = true;
                end
                if isfield(entity.(site_list(site).chain).(resfield),'CA')
                    CA_index = entity.(site_list(site).chain).(resfield).CA.tab_indices(c);
                else
                    failed = true;
                end
                if isfield(entity.(site_list(site).chain).(resfield_prev),'C')
                    C_index = entity.(site_list(site).chain).(resfield_prev).C.tab_indices(c);
                else
                    failed = true;
                end
                if failed % if not all aroms for proton construction were found
                    warnings = warnings + 1;
                    exceptions{warnings} = MException('pre_of_entity:no_proton',...
                        'For residue %i in chain %s proton is missing and cannot be constructed',site_list(site).residue,site_list(site).chain);
                    continue
                end
                % construct proton
                CN = entity.xyz(C_index,:)-entity.xyz(N_index,:);
                CAN = entity.xyz(CA_index,:)-entity.xyz(N_index,:);
                dir_HN = -(CN/norm(CN)+CAN/norm(CAN));
                dir_HN = dir_HN/norm(dir_HN);
                H_pred = entity.xyz(N_index,:) + NHbond*dir_HN;
                n_H = n_H + 1;
                coor(n_H,:) = H_pred;
                pre_list(n_H).chain = site_list(site).chain;
                pre_list(n_H).residue = site_list(site).residue;
            end
        end
    end
    
    % restrict lists and coordinates to number of actually found or constructed
    % protons
    pre_list = pre_list(1:n_H);
    coor = coor(1:n_H,:);
    
    label_adr = sprintf('{%i}%s',c,label_adr);
    [positions,entity,exceptions] = get_label(entity,label,'positions',label_adr);
    curr_label.coor = positions{1};
    [populations,entity] = get_label(entity,label,'populations',label_adr);
    curr_label.pop = populations{1};
    all_Gamma2 = pre(coor,curr_label,td,taur,taui,R2_dia,larmor);
    Gamma2(1:n_H) = Gamma2(1:n_H) + entity.populations(c)*all_Gamma2;
end

% copy the ensemble-average PREs to the output list
for k_H = 1:n_H
    pre_list(k_H).Gamma2 = Gamma2(k_H);
    pre_list(k_H).pre = R2_dia*exp(-td*Gamma2(k_H))/(R2_dia+Gamma2(k_H));
end