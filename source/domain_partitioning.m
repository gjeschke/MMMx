function entity = domain_partitioning(entity,options)
%
% DOMAIN_PARTITIONING   Partitions an AlphaFold entity into folded domains
%                       and flexible linkers
%
%   entity = domain_partitioning(entity)
%   Augments the entity, writes a rigid-body PDB file, and the frame of a
%   RigiFlex modelling script
%
%   entity = domain_partitioning(entity,options)
%   Allows for control of partitioning and of output
%
% INPUT
% entity        MMMx:atomic entity, must contain a predicted aligned error
%               matrix in field .pae, as AlphaFold models do; the current
%               version supports only single-chain entities
% options       .rbfile     name of output rigid-body file, defaults to
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
%
% OUTPUT
% entity       entity structure in MMMx:atomic representation 
%              the input entity is augmented by
%              .domains         (n,2) list of n folded domains specified by
%                               residue ranges, can be empty
%              .variability     (n,1) mean predicted aligned error (Å) 
%                               within each folded domain
%              .maxpae          (n,1) maxumum predicted aligned error (Å)
%                               within each folded domain
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
% Copyright(c) 2023: Gunnar Jeschke

if ~exist('options','var') || ~isfield(options,'rbfile')
    options.rbfile = sprintf('%s_rigid_bodies.pdb',entity.uniprot);
end

if ~isfield(options,'script')
    options.script = sprintf('%s_RigiFlex.mcx',entity.uniprot);
end

if ~isfield(options,'threshold') || isempty(options.threshold)
    options.threshold = 10;
end

if ~isfield(options,'unfify') || isempty(options.unify)
    options.unify = 15;
end

if ~isfield(options,'minsize') || isempty(options.minsize)
    options.minsize = 50;
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
    options.figname = sprintf('%s_domains.pdf',entity.uniprot);
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
    options.minrot = 5;
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

pairdist = (entity.pae + entity.pae')/2;

h = figure; clf; hold on
image(entity.pae,'CDataMapping','scaled');
curr_axis = gca;
set(curr_axis,'YDir','normal');
colorbar;
axis tight
xlabel('Residue number');
ylabel('Residue number');
title(sprintf('PAE and domains for %s',entity.uniprotname),'Interpreter','none');
axis equal

[n,~] = size(pairdist);

A = pairdist < options.threshold;

for k1 = 1:n-2*options.local
    for k2 = k1-options.local:k1+options.local
        if k2 >=1 && k2 <=n
            A(k1,k2) = 0;
            A(k2,k1) = 0;
        end
    end
end


% find diagonal blocks, based on F. Pedroche, M. Rebollo, C. Carrascosa and A. Palomares (2012)
% http://arxiv.org/abs/1206.5726
L = diag( sum(A,2) ) - A;
value = sum( triu(L) );
freq = find( value == 0 );


domains = zeros(length(freq),2);
variability = zeros(length(freq),1);
maxpae = zeros(length(freq),1);
dpoi = 0;
dstart = 1;
if length(freq) ==1
   mean_uncert = mean(mean(pairdist)); 
   if mean_uncert <= options.threshold
       dpoi = dpoi + 1;
       variability(dpoi) = mean_uncert;
       domains(dpoi,1) = dstart;
       domains(dpoi,2) = freq(1);
       maxpae(dpoi) =  max(max(pairdist(dstart:freq(1),dstart:freq(1))));
   end
else
    for k = 2:length(freq)
        if freq(k)-dstart >= options.minsize
            mean_uncert = mean(mean(pairdist(dstart:freq(k),dstart:freq(k))));
            if mean_uncert <= options.threshold
                dpoi = dpoi + 1;
                variability(dpoi) = mean_uncert;
                domains(dpoi,1) = dstart;
                domains(dpoi,2) = freq(k);
                maxpae(dpoi) =  max(max(pairdist(dstart:freq(k),dstart:freq(k))));
            end
        end
        dstart = freq(k);
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
                mean_uncert = mean(mean(pairdist(domains(cdpoi,1):domains(cdpoi,2),domains(cdpoi,1):domains(cdpoi,2))));
                variability(cdpoi) = mean_uncert;
                maxpae(cdpoi) = max(max(pairdist(domains(cdpoi,1):domains(cdpoi,2),domains(cdpoi,1):domains(cdpoi,2))));
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
    elseif max_terminal < options.minterminal
        entity.class = 'short flexible terminus';
    else
        entity.class = 'long flexible terminus';
    end
end


entity.domains = domains(1:dpoi,:);
entity.variability = variability(1:dpoi);
entity.maxpae = maxpae(1:dpoi);

if dpoi == 0
    lpoi = 1;
    linkers(1).range = [1,n];
    linkers(1).sequence = entity.sequence;
    domains = [];
else
    linkers(dpoi+1).range = [];
    linkers(dpoi+1).sequence = '';
    lpoi = 0;
    if domains(1,1) > 1
        lpoi = lpoi + 1;
        linkers(lpoi).range = [1,domains(1,1)-1];
        linkers(lpoi).sequence = entity.sequence(1:domains(1,1)-1);
    end
    for k = 2:dpoi
        lpoi = lpoi + 1;
        linkers(lpoi).range = [domains(k-1,2)+1,domains(k,1)-1];
        linkers(lpoi).sequence = entity.sequence(domains(k-1,2)+1:domains(k,1)-1);
    end
    if domains(dpoi,2) < n
        lpoi = lpoi + 1;
        linkers(lpoi).range = [domains(dpoi,2)+1,n];
        linkers(lpoi).sequence = entity.sequence(domains(dpoi,2)+1:n);
    end
end
entity.linkers = linkers(1:lpoi);

entity.folded = folded/n;

if~isempty(options.rbfile) && dpoi > 0
    mentity = rigid_domains(entity,domains(1:dpoi,:));
    put_pdb(mentity,options.rbfile);
end

if ~isempty(options.script)
    fid = fopen(options.script,'wt');
    fprintf(fid,'%% RigiFlex script template for AlphaFold prediction %s (%s)\n\n',entity.uniprot,entity.uniprotname);
    fprintf(fid,'#log\n\n');
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
                else
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
            fprintf(fid,'!flex %% N-terminal section\n');
            if rba
                fprintf(fid,'   expand %s_rba\n',entity.uniprot);
            else
                fprintf(fid,'   addpdb %s\n',options.rbfile);
            end
            modelfile = sprintf('%s_N_terminal',entity.uniprot);
            fprintf(fid,'   save %s\n',modelfile);
            fprintf(fid,'   sequence 1 %i %s\n',domains(1,1)-1,entity.sequence(1:domains(1,1)-1));
            fprintf(fid,'   c_anchor   (A)%i\n',domains(1,1));
            fprintf(fid,'.flex\n');
        end
        % make flex blocks for peptide linkers
        for k = 2:dpoi
            ctag1 = char(double('A')+k-2); % chain tag for first anchor domain
            ctag2 = char(double('A')+k-1); % chain tag for second anchor domain
            fprintf(fid,'\n');
            fprintf(fid,'!flex %% linker %i\n',k-1);
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
            fprintf(fid,'.flex\n');
        end
        % make flex block for C-terminal section if present and requested
        if domains(dpoi,2) < n && options.termini
            fprintf(fid,'\n');
            ctag = char(double('A')+dpoi-1);
            fprintf(fid,'!flex %% C-terminal section\n');
            if rba || N_terminal
                fprintf(fid,'   addpdb %s*.pdb\n',modelfile);
            else
                fprintf(fid,'   addpdb %s\n',options.rbfile);
            end
            modelfile = sprintf('%s_C_terminal',entity.uniprot);
            fprintf(fid,'   save %s\n',modelfile);
            fprintf(fid,'   sequence %i %i %s\n',domains(dpoi,2)+1,n,entity.sequence(domains(dpoi,2)+1:n));
            fprintf(fid,'   n_anchor   (%s)%i\n',ctag,domains(dpoi,2));
            fprintf(fid,'.flex\n');
        end
    end
    fprintf(fid,'\n#report\n');
    fclose(fid);
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