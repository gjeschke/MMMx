function entity = domain_partitioning(entity,options,eau)
%
% DOMAIN_PARTITIONING   Partitions an ensemble entity into folded domains
%                       and flexible linkers
%
%   entity = domain_partitioning(entity)
%   Augments the entity, writes a rigid-body PDB file, and (optionally) 
%   the frame of a RigiFlex modelling script
%
%   entity = domain_partitioning(entity,options,eau)
%   Allows for control of partitioning and of output
%
%   WARNING:    this is currently implemented only for entities consisting
%               of a single protein chain
%
% INPUT
% entity        MMMx:atomic entity; the current version supports only
%               single-chain entities
% options       .eau        flag. if true, an ensemble aligned uncertainty 
%                           matrix is generated and domain partitioning is 
%                           performed on this EAU matrix
%                           defaults to false if the entity contains an
%                           AlphaFold PAE matrix
%                           defaults to true, if an EAU matrix is provided
%                           as an argument, this overrides the AF default
%                           defaults to true, if the entity does not
%                           contain a PAE matrix and has several
%                           conformers
%               .rbfile     name of output rigid-body file, defaults to
%                           <entity.name>_rigid_bodies.pdb, if present but
%                           empty, the rigid-body file is not written
%               .script     name of output RigiFlex script, defaults to
%                           <entity.name>_RigiFlex.mcx, if present but
%                           empty, no RigiFlex script is written
%               .threshold  maximum predicted aligned error inside domain,
%                           defaults to 10 Angstroem
%               .unify      maximum predicted aligned error for unifying
%                           domains connected by short linkers, defaults to
%                           15 Angstroem
%               .minsize    minimum domain size, defaults to 50 residues
%               .minlink    minimum linker length, defaults to 3 residues
%               .minterm    minimum terminal section length, used only for
%                           classification, defaults to 10
%               .local      local context length, defaults to 5 residues
%               .figname    name of output figure file, defaults to
%                           <entity.name>_domains.pdf, if present and
%                           empty, figure is only displayed, not saved
%               .rmin       minimum distance for rigid-body reference
%                           points, defaults to 15 Å
%               .rmax       maximum distance for rigid-body reference
%                           points, defaults to 80 Å
%               .restypes   residue types allowed for spin labelling,
%                           defaults to 'CILMSTVAN'
%               .minrot     minimum number of rotamers in labelling,
%                           defaults to 5
%               .minZ       minimum partition function in labelling,
%                           defaults to 0.1
%               .label      spin label, defaults to 'mtsl'
%               .termini    generate flex blocks for terminal sections,
%                           defaults to true
% eau           ensemble aligned uncertainty matrix
%
% OUTPUT
% entity       entity structure in MMMx:atomic representation 
%              the input entity is augmented by
%              .domains         (n,2) list of n folded domains specified by
%                               residue ranges, can be empty
%              .linkers         list of m flexible linkers, each with
%                               subfields
%                               .range      (1,2) first, second residue
%                               .sequence   amino acid sequence
%              .folded          fraction of residues in folded domains
%              .class           flexibility class, one of the classes
%                               'folded'
%                               'short flexible terminus'
%                               'long flexible terminus'
%                               'RigiFlex'
%                               'disordered'
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023-2025: Gunnar Jeschke

cap = 31.75; % cap of AlphaFold2 PAE

PAE2stddev = 1; % factor for converting PAE to standard deviation 

my_base = load('deer_kernel.mat','base','t','r');
Pake_base.kernel = my_base.base-ones(size(my_base.base)); % kernel format for pcf2deer
Pake_base.t = my_base.t;
Pake_base.r = my_base.r;
clear my_base

if ~exist('options','var') || ~isfield(options,'rbfile')
    options.rbfile = sprintf('%s_rigid_bodies.pdb',entity.uniprot);
end

if ~isfield(options,'script')
    options.script = sprintf('%s_RigiFlex.mcx',entity.uniprot);
    options.fit_script = sprintf('%s_EnsembleFit.mcx',entity.uniprot);
elseif ~isempty(options.script)
    [pname,bname] = fileparts(options.script);
    options.fit_script = fullfile(pname,sprintf('%s_EnsembleFit.mcx',bname));
end

if ~isfield(options,'threshold') || isempty(options.threshold)
    options.threshold = cap/3;
end

if ~isfield(options,'unfify') || isempty(options.unify)
    options.unify = 15;
end

if ~isfield(options,'minsize') || isempty(options.minsize)
    options.minsize = 25;
end

if ~isfield(options,'minlink') || isempty(options.minlink)
    options.minlink = 3;
end

if ~isfield(options,'local') || isempty(options.local)
    options.local = 5;
end

if ~isfield(options,'minterm') || isempty(options.minterm)
    options.minterminal = 10;
end

if ~isfield(options,'figname')
    options.figname = sprintf('%s_domains.pdf',entity.uniprotname);
end

if ~isfield(options,'rmin') || isempty(options.rmin)
    options.rmin = 15;
end

if ~isfield(options,'rmax') || isempty(options.rmax)
    options.rmax = 80;
end

if ~isfield(options,'restypes') || isempty(options.restypes)
    options.restypes = 'CILMSTVAN';
end

if ~isfield(options,'minrot') || isempty(options.minrot)
    if isfield(options,'label') && length(options.label) > 5 && strcmpi(options.label(1:5),'atom.')
        options.minrot = 1;
    else
        options.minrot = 5;
    end
end

if ~isfield(options,'minZ') || isempty(options.minZ)
    options.minZ = 0.1;
end

if ~isfield(options,'label') || isempty(options.label)
    options.label = 'mtsl';
end

if ~isfield(options,'termini') || isempty(options.termini)
    options.termini = true;
end

if ~isfield(options,'eau')
    if isfield(entity,'pae')
        options.eau = false;
    elseif length(entity.populations) > 1
        options.eau = true;
    end
    if exist('eau','var')
        options.eau = true;
    end
end

h = figure; clf; hold on

if options.eau
    if ~isfield(entity,'eau')
        if exist('eau','var')
            entity.eau = eau;
        else
            entity = ensemble_aligned_uncertainty(entity);
        end
    end
    pairdist = (entity.eau + entity.eau')/2;
    atype = 'EAU';
    image(entity.eau,'CDataMapping','scaled');
else
    pairdist = (entity.pae + entity.pae')/2;
    atype = 'PAE';
    image(entity.pae,'CDataMapping','scaled');
end

% analyze the first chain
chains = fieldnames(entity);
for c = 1:length(chains)
    chain = chains{c};
    if isstrprop(chain(1),'upper')
        break
    end
end

curr_axis = gca;
set(curr_axis,'YDir','normal');
colorbar;
axis tight
xlabel('Residue number');
ylabel('Residue number');
title(sprintf('%s and domains for %s',atype,entity.name),'Interpreter','none');
axis equal

[n,~] = size(pairdist);

A = pairdist < options.threshold;

extension = ones(1,n);
for k1 = 2:n-1
    ext = 1;
    k2 = k1+options.local;
    while k2 < n
        k2 = k2 + 1;
        if A(k1,k2)
            ext = k2-k1;
        else
            break
        end
    end
    extension(k1) = ext;
end
% moving average filter
b = ones(1,options.local)/options.local;
a = 1;
extension = filter(b,a,extension);

dpoi = 0;
domains = zeros(ceil(n/options.minsize),2);
while max(extension) >= options.minsize
    [ext,k] = max(extension);
    dstart = k;
    dend = k + round(ext);
    if dend > n
        dend = n;
    end
    mean_uncert0 = mean(mean(pairdist(dstart:dend,dstart:dend)));
    for k1 = dstart:-1:1
        mean_uncert = mean(pairdist(k1,k1:dend));
        if mean_uncert <= options.threshold
            dstart = k1;
            mean_uncert0 = mean(mean(pairdist(dstart:dend,dstart:dend)));
        else
            break;
        end
    end
    if mean_uncert0 > options.threshold
        for k1 = dstart:n
            mean_uncert = mean(mean(pairdist(k1:dend,k1:dend)));
            if mean_uncert < mean_uncert0
                mean_uncert0 = mean_uncert;
                dstart = k1;
            else
                break;
            end
        end
    end
    for k1 = dend:n
        mean_uncert = mean(mean(pairdist(k1,dstart:k1)));
        if mean_uncert <=options.threshold
            mean_uncert0 = mean_uncert;
            dend = k1;
        else
            break;
        end
    end
    if mean_uncert0 > options.threshold
        for k1 = dend:-1:1
            mean_uncert = mean(mean(pairdist(dstart:k1,dstart:k1)));
            if mean_uncert < mean_uncert0
                mean_uncert0 = mean_uncert;
                dend = k1;
            else
                break;
            end
        end
    end
    extension(dstart:dend) = zeros(1,dend-dstart+1);
    pairdist(dstart:dend,dstart:dend) = 1e12;
    for k = 1:n
        test = k + round(extension(k));
        if test <= n && extension(test) > cap
            extension(k) = 1;
        end
    end
    if mean_uncert0 <= options.threshold
        dpoi = dpoi +1;
        domains(dpoi,:) = [dstart,dend];
    end
end

% unify domains separated by too short linkers, if requested
if options.minlink > 1
    cdpoi = 1;
    for k = 2:dpoi
        if domains(k,1) - domains(k-1,2) < options.minlink + 1
            test = max(max(pairdist(domains(k-1,1):domains(k,2),domains(k-1,1):domains(k,2))));
            if test <= options.unify
                domains(cdpoi,2) = domains(k,2);
            elseif domains(k,1) - domains(k-1,2) < 3
                domains(k,1) = domains(k,1) + 1;
                domains(k-1,2) = domains(k-1,2) - 1;
                cdpoi = cdpoi + 1;
            else
                cdpoi = cdpoi + 1;
            end
        else
            cdpoi = cdpoi + 1;
        end
    end
    dpoi = cdpoi;
end

if max(max(domains)) == 0
    dpoi = 0;
else
    % sort domains
    [~,indices] = sort(domains(1:dpoi,1));
    domains = domains(indices,:);
end


folded = 0;
for k = 1:dpoi
    k1 = domains(k,1);
    k2 = domains(k,2);
    folded = folded + k2 - k1 + 1;
    % fprintf(1,'Folded domain: (%i,%i)\n',k1,k2);
    plot([k1,k1],[k1,k2],'LineWidth',2,'Color',[0.8,0,0]);
    plot([k2,k2],[k1,k2],'LineWidth',2,'Color',[0.8,0,0]);
    plot([k1,k2],[k1,k1],'LineWidth',2,'Color',[0.8,0,0]);
    plot([k1,k2],[k2,k2],'LineWidth',2,'Color',[0.8,0,0]);
end

if ~isempty(options.figname)
    saveas(h,options.figname);
end

max_terminal = NaN;
if dpoi > 0
    max_terminal = domains(1,1) - 1;
    if n - domains(dpoi,2) > max_terminal
        max_terminal = n - domains(end,2);
    end
end

% classify this protein
if dpoi < 1
    entity.class = 'disordered';
else
    if dpoi > 1
        entity.class = 'RigiFlex';
    elseif max_terminal == 0
        entity.class = 'folded';
    elseif max_terminal < options.minterm
        entity.class = 'short flexible terminus';
    else
        entity.class = 'long flexible terminus';
    end
end


entity.domains = domains(1:dpoi,:);

if ~isfield(entity.(chain),'sequence')
    entity = add_sequence_information(entity,chain);
end

if dpoi == 0
    lpoi = 1;
    linkers(1).range = [1,n];
    linkers(1).sequence = entity.(chain).sequence;
    domains = [];
else
    linkers(dpoi+1).range = [];
    linkers(dpoi+1).sequence = '';
    lpoi = 0;
    if domains(1,1) > 1
        lpoi = lpoi + 1;
        linkers(lpoi).range = [1,domains(1,1)-1];
        linkers(lpoi).sequence = entity.(chain).sequence(1:domains(1,1)-1);
    end
    for k = 2:dpoi
        lpoi = lpoi + 1;
        linkers(lpoi).range = [domains(k-1,2)+1,domains(k,1)-1];
        linkers(lpoi).sequence = entity.(chain).sequence(domains(k-1,2)+1:domains(k,1)-1);
    end
    if domains(dpoi,2) < n
        lpoi = lpoi + 1;
        linkers(lpoi).range = [domains(dpoi,2)+1,n];
        linkers(lpoi).sequence = entity.(chain).sequence(domains(dpoi,2)+1:n);
    end
end
entity.linkers = linkers(1:lpoi);

entity.folded = folded/n;

if~isempty(options.rbfile) && dpoi > 0
    mentity = rigid_domains(entity,domains(1:dpoi,:));
    put_pdb(mentity,options.rbfile);
    if strcmpi(options.label,'atom.CA')
        mentity = get_pdb(options.rbfile);
    end
end

if ~isempty(options.script)
    fid = fopen(options.script,'wt');
    fprintf(fid,'%% RigiFlex script template for AlphaFold prediction %s (%s)\n\n',entity.uniprot,entity.uniprotname);
    fprintf(fid,'#log\n\n');
    fidef = fopen(options.fit_script,'wt');
    fprintf(fidef,'%% EnsembleFit script template for AlphaFold prediction %s (%s)\n\n',entity.uniprot,entity.uniprotname);
    fprintf(fidef,'#log\n\n');
    fprintf(fidef,'!EnsembleFit\n');
    fprintf(fidef,'   csv\n');
    fprintf(fidef,'   plot\n');
    fprintf(fidef,'   nnllsq\n');
    fprintf(fidef,'   save %s_ensemble.ens',entity.uniprot);
    fprintf(fidef,'\n   deer atom.CA\n');
    rba = false;
    if dpoi > 1 % Rigi block is written only if more than one folded domain exists
        rba = true;
        fprintf(fid,'!rigi\n');
        fprintf(fid,'   rbtemplate %s\n',options.rbfile);
        fprintf(fid,'   separate on\n');
        fprintf(fid,'   maxtrials 10000\n');
        fprintf(fid,'   models 200\n');
        fprintf(fid,'   save %s_rba\n',entity.uniprot);
        % generate reference points for rigid bodies
        rb_entity = get_pdb(options.rbfile);
        chains = fieldnames(rb_entity);
        refpoint_options.r_min = options.rmin;
        refpoint_options.r_max = options.rmax;
        sitescan_options.restypes = options.restypes;
        sitescan_options.min_rotamers = options.minrot;
        sitescan_options.min_Z = options.minZ;
        sitescan_options.ensemble = false;
        rbi = 0;
        all_refpoints = cell(length(chains),1);
        complete_reference = true;
        for kc = 1:length(chains)
            chain = chains{kc};
            if isstrprop(chain(1),'upper') % chain fields start with a capital
                rbi = rbi + 1;
                sitescan_options.chains = chain;
                ssfname = sprintf('%s_%s_sites',rb_entity.name,chain);
                rb_entity = site_scan_labels(rb_entity,options.label,ssfname,sitescan_options);
                all_site_coor_pop(10000).coor = [];
                all_site_coor_pop(10000).pop = [];
                all_site_coor_pop(10000).label = '';
                all_site_coor_pop(10000).address = '';
                [sites,labels] = rd_site_list(ssfname);
                [rb_entity,site_coor_pop] = site_geometries(rb_entity,sites,labels);
                [refpoints,failed] = get_refpoints(site_coor_pop,refpoint_options);
                if failed
                    fprintf(fid,'\nReference point computation failed for rigid body %s\n',chain);
                    complete_reference = false;
                else
                    all_refpoints{kc} = refpoints;
                    fprintf(fid,'   rigid (%s) %% %i-%i\n',chain,domains(rbi,:));
                    for k = 1:3
                        fprintf(fid,'      %s %s\n',refpoints(k).site,refpoints(k).label);
                    end
                    fprintf(fid,'   .rigid\n\n');
                end
            end
        end
        % make linker statements
        for k = 2:dpoi
            ctag1 = char(double('A')+k-2); % chain tag for first anchor domain
            ctag2 = char(double('A')+k-1); % chain tag for second anchor domain
            fprintf(fid,'\n');
            fprintf(fid,'   plink %% linker %i\n',k-1);
            fprintf(fid,'      (%s)%i   (%s)%i   %i\n',ctag1,domains(k-1,2),ctag2,domains(k,1),domains(k,1)-domains(k-1,2));
            fprintf(fid,'   .plink\n');
        end
        % make distance distribution restraints, if AlphaFold-based modeling is possible and requested 
        if complete_reference && strcmpi(options.label,'atom.CA')
            site_pair_list = zeros(100,2);
            spl_pointer = 0;
            fprintf(fid,'\n   ddr atom.CA\n');
            % restraints between reference points
            nr1 = 0;
            for kc1 = 1:length(chains)-1
                chain = chains{kc1};
                if isstrprop(chain(1),'upper') % chain fields start with a capital
                    nr1 = nr1 + 1;
                    ref1 = all_refpoints{kc1};
                    for kr1 = 1:3
                        coor1 = get_label(mentity,ref1(kr1).label,'coor',ref1(kr1).site);
                        poi = strfind(ref1(kr1).site,')');
                        res1 = str2double(ref1(kr1).site(poi+1:end));
                        poi = strfind(ref1(kr1).site,'}');
                        if isempty(poi)
                            poi = 0;
                        end
                        site1 = ref1(kr1).site(poi+1:end);
                        nr2 = nr1;
                        for kc2 = kc1+1:length(chains)
                            chain = chains{kc2};
                            if isstrprop(chain(1),'upper') % chain fields start with a capital
                                nr2 = nr2 + 1;
                                ref2 = all_refpoints{kc2};
                                for kr2 = 1:3
                                    coor2 = get_label(mentity,ref2(kr2).label,'coor',ref2(kr2).site);
                                    poi = strfind(ref2(kr2).site,')');
                                    res2 = str2double(ref2(kr2).site(poi+1:end));
                                    poi = strfind(ref2(kr2).site,'}');
                                    if isempty(poi)
                                        poi = 0;
                                    end
                                    site2 = ref2(kr2).site(poi+1:end);
                                    rmean = norm(coor1{1}-coor2{1});
                                    sigr = max([entity.pae(res1,res2),entity.pae(res2,res1)])*PAE2stddev; 
                                    if sigr <= 30*PAE2stddev
                                        fprintf(fid,'      %s   %s   %4.1f   %4.1f %% ref %i.%i to %i.%i\n',site1,site2,rmean,sigr,nr1,kr1,nr2,kr2);
                                        sim_deer_file = sprintf('sim_deer_CA_CA_%s_%s.csv',site1,site2);
                                        [t,ff] = get_ff(rmean,sigr,Pake_base);
                                        fprintf(fidef,'      (A)%i   (A)%i   @%s %% ref %i.%i to %i.%i\n',res1,res2,sim_deer_file,nr1,kr1,nr2,kr2);
                                        sim_deer = [t',ff',ff',zeros(size(ff'))];
                                        description = {'% time (ns)','simulated','simulated','zero background'};
                                        put_csv(sim_deer_file,sim_deer,description);
                                        spl_pointer = spl_pointer + 1;
                                        site_pair_list(spl_pointer,:) = [res1,res2];
                                    end
                                end
                            end
                        end
                    end
                end
            end
            % auxiliary restraints
            for d1 = 1:dpoi-1
                d1_offset = domains(d1,1)-1;
                chain1 = char(double('A') -1 + d1);
                for d2 = d1+1:dpoi
                    chain2 = char(double('A') -1 + d2);
                    d2_offset = domains(d2,1)-1;
                    id_pae_upper = entity.pae(domains(d1,1):domains(d1,2),domains(d2,1):domains(d2,2));
                    id_pae_lower = entity.pae(domains(d2,1):domains(d2,2),domains(d1,1):domains(d1,2))';
                    id_sigr = id_pae_upper;
                    id_sigr(id_sigr < id_pae_lower) = id_pae_lower(id_sigr < id_pae_lower);
                    id_sigr = id_sigr*PAE2stddev; 
                    min_err = min(min(id_sigr));
                    aux = 0;
                    while min_err <= 30*PAE2stddev && aux < 20
                        [colmin,rows] = min(id_sigr);
                        [min_err,col] = min(colmin);
                        row = rows(col);
                        add_it = true;
                        for sp = 1:spl_pointer % search for existing sitepairs very close to this one
                            if abs(site_pair_list(sp,1) - (row + d1_offset)) < 5 &&  abs(site_pair_list(sp,2) - (col + d2_offset)) < 5
                                add_it = false;
                            end
                        end
                        id_sigr(row,col) = 1000;
                        coor1 =  get_label(entity,'atom.CA','coor',sprintf('(A)%i',row + d1_offset));
                        coor2 =  get_label(entity,'atom.CA','coor',sprintf('(A)%i',col + d2_offset));
                        rmean = norm(coor1{1}-coor2{1});
                        if add_it
                            aux = aux + 1;
                            spl_pointer = spl_pointer + 1;
                            site_pair_list(spl_pointer,:) = [row + d1_offset,col + d2_offset];
                            fprintf(fid,'      (%c)%i   (%c)%i   %4.1f   %4.1f %% aux %i to %i\n',chain1,row + d1_offset,chain2,col + d2_offset,rmean,min_err,d1,d2);
                            sim_deer_file = sprintf('sim_deer_CA_CA_(%c)%i_(%c)%i.csv',chain1,row + d1_offset,chain2,col + d2_offset);
                            [t,ff] = get_ff(rmean,sigr,Pake_base);
                            fprintf(fidef,'      (A)%i   (A)%i   @%s %% aux %i to %i\n',row + d1_offset,col + d2_offset,sim_deer_file,d1,d2);
                            sim_deer = [t',ff',ff',zeros(size(ff'))];
                            description = {'% time (ns)','simulated','simulated','zero background'};
                            put_csv(sim_deer_file,sim_deer,description);
                        end
                        min_err = min(min(id_sigr));
                    end
                end
            end
            fprintf(fid,'   .ddr\n');
        end
        fprintf(fid,'\n.rigi\n');
    end
    if dpoi == 0
        fprintf(fid,'\n');
        fprintf(fid,'!flex %% whole protein\n');
        modelfile = sprintf('%s_full_length',entity.uniprot);
        fprintf(fid,'   save %s\n',modelfile);
        fprintf(fid,'   sequence 1 %i %s\n',n,entity.sequence(1:n));
        fprintf(fid,'.flex\n');
    else
        N_terminal = false;
        % make flex block for N-terminal section if present and requested
        if domains(1,1) > 1 && options.termini
            N_terminal = true;
            fprintf(fid,'\n');
            fprintf(fid,'!flex 0.5 1 0.25 %% N-terminal section\n');
            if rba
                fprintf(fid,'   expand %s_rba\n',entity.uniprot);
            else
                fprintf(fid,'   addpdb %s\n',options.rbfile);
            end
            modelfile = sprintf('%s_N_terminal',entity.uniprot);
            fprintf(fid,'   save %s\n',modelfile);
            fprintf(fid,'   sequence 1 %i %s\n',domains(1,1)-1,entity.sequence(1:domains(1,1)-1));
            fprintf(fid,'   c_anchor   (A)%i\n',domains(1,1));
            % make distance distribution restraints, if requested
            if strcmpi(options.label,'atom.CA')
                fprintf(fid,'\n   ddr atom.CA\n');
                for refd = 1:dpoi % loop over reference domains
                    chain1 = char(double('A') -1 + refd);
                    unrestrained = 0;
                    for res = 1:domains(1,1)-1 % loop over residues in the N-terminal linker
                        id_pae_upper = entity.pae(res,domains(refd,1):domains(refd,2));
                        id_pae_lower = entity.pae(domains(refd,1):domains(refd,2),res)';
                        id_sigr = id_pae_upper;
                        id_sigr(id_sigr < id_pae_lower) = id_pae_lower(id_sigr < id_pae_lower);
                        id_sigr = id_sigr*PAE2stddev; 
                        [min_err,refpoi] = min(id_sigr);
                        if min_err <= 30*PAE2stddev
                            if min_err/unrestrained <= 10*PAE2stddev
                                coor1 =  get_label(entity,'atom.CA','coor',sprintf('(A)%i',refpoi + domains(refd,1)-1));
                                coor2 =  get_label(entity,'atom.CA','coor',sprintf('(A)%i',res));
                                rmean = norm(coor1{1}-coor2{1});
                                fprintf(fid,'      (%c)%i   %i   %4.1f   %4.1f %% N-terminal flexible domain to folded domain %i\n',chain1,refpoi + domains(refd,1)-1,res,rmean,min_err,refd);
                                sim_deer_file = sprintf('sim_deer_CA_CA_(A)%i_(A)%i.csv',refpoi + domains(refd,1)-1,res);
                                [t,ff] = get_ff(rmean,min_err,Pake_base);
                                sim_deer = [t',ff',ff',zeros(size(ff'))];
                                description = {'% time (ns)','simulated','simulated','zero background'};
                                put_csv(sim_deer_file,sim_deer,description);
                                fprintf(fidef,'      (A)%i   (A)%i   @%s %% N-terminal loop to domain %i\n',refpoi + domains(refd,1)-1,res,sim_deer_file,refd);
                                unrestrained = 0;
                            else
                                unrestrained = unrestrained + 1;
                            end
                        else
                            unrestrained = unrestrained + 1;
                        end
                    end
                end
                fprintf(fid,'   .ddr\n');
            end
            fprintf(fid,'.flex\n');
        end
        % make flex blocks for peptide linkers
        for k = 2:dpoi
            ctag1 = 'A'; % chain tag for first anchor domain
            ctag2 = char(double('A')+k-1); % chain tag for second anchor domain
            fprintf(fid,'\n');
            fprintf(fid,'!flex 0.5 1 0.25 %% linker %i\n',k-1);
            if k == 2 && ~N_terminal
                fprintf(fid,'   expand %s_rba\n',entity.uniprot);
            else
                fprintf(fid,'   addpdb %s*.pdb\n',modelfile);
            end
            modelfile = sprintf('%s_linker_%i',entity.uniprot,k-1);
            fprintf(fid,'   save %s\n',modelfile);
            fprintf(fid,'   sequence %i %i %s\n',domains(k-1,2)+1,domains(k,1)-1,entity.sequence(domains(k-1,2)+1:domains(k,1)-1));
            fprintf(fid,'   n_anchor   (%s)%i\n',ctag1,domains(k-1,2));
            fprintf(fid,'   c_anchor   (%s)%i\n',ctag2,domains(k,1));
            % make distance distribution restraints, if requested
            if strcmpi(options.label,'atom.CA')
                fprintf(fid,'\n   ddr atom.CA\n');
                for refd = 1:dpoi % loop over reference domains
                    if refd < k
                        chain1 = 'A';
                    else
                        chain1 = char(double('A') + refd - 1);
                    end
                    unrestrained = 0;
                    for res = domains(k-1,2)+1:domains(k,1)-1 % loop over residues in this flexible linker
                        id_pae_upper = entity.pae(res,domains(refd,1):domains(refd,2));
                        id_pae_lower = entity.pae(domains(refd,1):domains(refd,2),res)';
                        id_sigr = id_pae_upper;
                        id_sigr(id_sigr < id_pae_lower) = id_pae_lower(id_sigr < id_pae_lower);
                        id_sigr = id_sigr*PAE2stddev; 
                        [min_err,refpoi] = min(id_sigr);
                        if min_err <= 30*PAE2stddev
                            if min_err/unrestrained <= 10*PAE2stddev
                                coor1 =  get_label(entity,'atom.CA','coor',sprintf('(A)%i',refpoi + domains(refd,1)-1));
                                coor2 =  get_label(entity,'atom.CA','coor',sprintf('(A)%i',res));
                                rmean = norm(coor1{1}-coor2{1});
                                fprintf(fid,'      (%c)%i   %i   %4.1f   %4.1f %% linker %i to domain %i\n',chain1,refpoi + domains(refd,1)-1,res,rmean,min_err,k-1,refd);
                                sim_deer_file = sprintf('sim_deer_CA_CA_(A)%i_(A)%i.csv',refpoi + domains(refd,1)-1,res);
                                [t,ff] = get_ff(rmean,min_err,Pake_base);
                                sim_deer = [t',ff',ff',zeros(size(ff'))];
                                description = {'% time (ns)','simulated','simulated','zero background'};
                                put_csv(sim_deer_file,sim_deer,description);
                                fprintf(fidef,'      (A)%i   (A)%i   @%s %% linker %i-%i to domain %i\n',refpoi + domains(refd,1)-1,res,sim_deer_file,k-1,k,refd);
                                unrestrained = 0;
                            else
                                unrestrained = unrestrained + 1;
                            end
                        else
                            unrestrained = unrestrained + 1;
                        end
                    end
                end
                fprintf(fid,'   .ddr\n');
            end
            fprintf(fid,'.flex\n');
        end
        % make flex block for C-terminal section if present and requested
        if domains(dpoi,2) < n && options.termini
            fprintf(fid,'\n');
            ctag = 'A';
            fprintf(fid,'!flex 0.5 1 0.25 %% C-terminal section\n');
            if rba || N_terminal
                fprintf(fid,'   addpdb %s*.pdb\n',modelfile);
            else
                fprintf(fid,'   addpdb %s\n',options.rbfile);
            end
            modelfile = sprintf('%s_C_terminal',entity.uniprot);
            fprintf(fid,'   save %s\n',modelfile);
            fprintf(fid,'   sequence %i %i %s\n',domains(dpoi,2)+1,n,entity.sequence(domains(dpoi,2)+1:n));
            fprintf(fid,'   n_anchor   (%s)%i\n',ctag,domains(dpoi,2));
            % make distance distribution restraints, if requested
            if strcmpi(options.label,'atom.CA')
                fprintf(fid,'\n   ddr atom.CA\n');
                for refd = 1:dpoi % loop over reference domains
                    chain1 = 'A';
                    unrestrained = 0;
                    for res = 1:domains(1,1)-1 % loop over residues in the N-terminal linker
                        id_pae_upper = entity.pae(res,domains(refd,1):domains(refd,2));
                        id_pae_lower = entity.pae(domains(refd,1):domains(refd,2),res)';
                        id_sigr = id_pae_upper;
                        id_sigr(id_sigr < id_pae_lower) = id_pae_lower(id_sigr < id_pae_lower);
                        id_sigr = id_sigr*PAE2stddev; 
                        [min_err,refpoi] = min(id_sigr);
                        if min_err <= 30*PAE2stddev
                            if min_err/unrestrained <= 10*PAE2stddev
                                coor1 =  get_label(entity,'atom.CA','coor',sprintf('(A)%i',refpoi + domains(refd,1)-1));
                                coor2 =  get_label(entity,'atom.CA','coor',sprintf('(A)%i',res));
                                rmean = norm(coor1{1}-coor2{1});
                                fprintf(fid,'      (%c)%i   %i   %4.1f   %4.1f %% C-terminal flexible domain to folded domain %i\n',chain1,refpoi + domains(refd,1)-1,res,rmean,min_err,refd);
                                sim_deer_file = sprintf('sim_deer_CA_CA_(A)%i_(A)%i.csv',refpoi + domains(refd,1)-1,res);
                                [t,ff] = get_ff(rmean,min_err,Pake_base);
                                sim_deer = [t',ff',ff',zeros(size(ff'))];
                                description = {'% time (ns)','simulated','simulated','zero background'};
                                put_csv(sim_deer_file,sim_deer,description);
                                fprintf(fidef,'      (A)%i   (A)%i   @%s %% C-terminal loop to domain %i\n',refpoi + domains(refd,1)-1,res,sim_deer_file,refd);
                                unrestrained = 0;
                            else
                                unrestrained = unrestrained + 1;
                            end
                        else
                            unrestrained = unrestrained + 1;
                        end
                    end
                end
                fprintf(fid,'   .ddr\n');
            end
            fprintf(fid,'.flex\n');
        end
    end
    fprintf(fid,'\n#report\n');
    fclose(fid);
    fprintf(fidef,'   .deer\n\n');
    if exist('modelfile','var')
        fprintf(fidef,'   addpdb %s*.pdb\n',modelfile);
    end
    fprintf(fidef,'.EnsembleFit\n');
    fprintf(fidef,'\n#report\n');
    fclose(fidef);
end

function mentity = rigid_domains(entity,domains)

[d,~] = size(domains); 
mentity.name = entity.name;
mentity.populations = 1;
mentity.selected = 1;
mentity.index_array = zeros(100000,5,'uint16');
xyz = zeros(100000,3);
elements = zeros(1,100000,'uint8');
occupancies = ones(100000,1,'uint8');
atnum = 0;
for p = 1:d
    ctag = char(double('A')+p-1); % chain tag for this part in new entity
    mentity.(ctag).selected = false;
    mentity.(ctag).index = p;
    conf = 1;
    chain = 'A';
    range = domains(p,:);
    for r = range(1):range(2) % loop over all selected residues
        residue = sprintf('R%i',r);
        mentity.(ctag).(residue).index = 1 + r - range(1);
        mentity.(ctag).(residue).selected = false;
        mentity.(ctag).(residue).selected_rotamers = 1;
        mentity.(ctag).(residue).name = entity.(chain).(residue).name;
        mentity.(ctag).(residue).locations = ' ';
        mentity.(ctag).(residue).populations = 1;
        atoms = fieldnames(entity.(chain).(residue));
        for ka = 1:length(atoms)
            atom = atoms{ka};
            if isstrprop(atom(1),'upper') % these are atom fields
                this_atom = entity.(chain).(residue).(atom);
                at_index = this_atom.tab_indices(conf);
                atnum = atnum + 1;
                this_atom.tab_indices = atnum;
                xyz(atnum,:) = entity.xyz(at_index,:);
                elements(atnum) = entity.elements(at_index);
                mentity.(ctag).(residue).(atom) = this_atom;
            end
        end
    end
end
mentity.xyz = xyz(1:atnum,:);
mentity.elements = elements(1:atnum);
mentity.occupancies = occupancies(1:atnum);
mentity.index_array = mentity.index_array(1:atnum,:);

function [refpoints,failed] = get_refpoints(all_site_coor_pop,refpoint_options)

S = length(all_site_coor_pop); % total number of sites
max_area = 0;
best_sites = zeros(1,3);
% check all triangles
for s1 = 1:S-2
    coor1 = all_site_coor_pop(s1).coor;
    pop1 = all_site_coor_pop(s1).pop;
    for s2 = s1+1:S-1
        coor2 = all_site_coor_pop(s2).coor;
        pop2 = all_site_coor_pop(s2).pop;
        a = get_r_mean(coor1,pop1,coor2,pop2);
        if a < refpoint_options.r_min || a > refpoint_options.r_max
            continue
        end
        for s3 = s2+1:S
            coor3 = all_site_coor_pop(s3).coor;
            pop3 = all_site_coor_pop(s3).pop;
            b = get_r_mean(coor1,pop1,coor3,pop3);
            if b < refpoint_options.r_min || b > refpoint_options.r_max
                continue
            end
            c = get_r_mean(coor2,pop2,coor3,pop3);
            if c < refpoint_options.r_min || c > refpoint_options.r_max
                continue
            end
            s = (a+b+c)/2; % half the circumference
            area = sqrt(s*(s-a)*(s-b)*(s-c)); % area of the triangle, Heron's theorem
            if area > max_area
                max_area = area;
                best_sites = [s1,s2,s3];
            end
        end
    end
end

refpoints(1).site = '';
refpoints(2).site = '';
refpoints(3).site = '';
refpoints(1).label = '';
refpoints(2).label = '';
refpoints(3).label = '';
if prod(best_sites) == 0
    failed = true;
    return
else
    failed = false;
    for k = 1:3
        refpoints(k).site = all_site_coor_pop(best_sites(k)).address;
        refpoints(k).label = all_site_coor_pop(best_sites(k)).label;
    end
end

function [c_entity,site_coor_pop] = site_geometries(c_entity,sites,label)

site_coor_pop(length(sites)).coor = [];
site_coor_pop(length(sites)).pop = [];
site_coor_pop(length(sites)).label = '';
site_coor_pop(length(sites)).address = '';

poi = 0;
for k = 1:length(sites)
    [argsout,c_entity,exceptions] = get_label(c_entity,label,{'positions','populations'},sites(k).address);
    if ~isempty(exceptions{1})
        continue
    end
    poi = poi + 1;
    site_coor_pop(poi).coor = argsout{1}{1};
    site_coor_pop(poi).pop = argsout{2}{1};
    site_coor_pop(poi).label = label;
    site_coor_pop(poi).address = sites(k).address;
end
site_coor_pop = site_coor_pop(1:poi);

function r_mean = get_r_mean(positions1,populations1,positions2,populations2)
% Computes a pair distance distribution as well as the fraction of missing
% population, the distribution is normalized to unity

% normalize population vectors
populations1 = populations1/sum(populations1);
populations2 = populations2/sum(populations2);

[n_rot_1,~] = size(positions1);   % returns number of rotamers at position 1
[n_rot_2,~] = size(positions2);   % -//- at position 2

a2 = repmat(sum(positions1.^2,2),1,n_rot_2);
b2 = repmat(sum(positions2.^2,2),1,n_rot_1).';
pair_dist = sqrt(abs(a2 + b2 - 2*positions1*positions2.'));
weights = populations1*populations2.';
r_mean = sum(sum(pair_dist.*weights));

function [t,sim_ff] = get_ff(rmean,sigr,Pake_base)

arg = 10*(Pake_base.r - rmean/10)/sigr;
distr = exp(-arg.^2/2);
distr = distr/sum(distr);
sim_ff = form_factor_deer(Pake_base.r,distr,Pake_base.t,Pake_base);
determine = abs(sim_ff) < 0.02;
n = length(determine);
while Pake_base.t(n) > 5 && determine(n) == 1
    n = n - 1;
end
sim_ff = sim_ff(1:n);
t = Pake_base.t(1:n);

