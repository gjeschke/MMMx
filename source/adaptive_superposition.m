function entity = adaptive_superposition(entity,options,eau)
%
% ADAPTIVE_SUPERPOSITION Tries to identify folded domains and superimposes
%             the ensemble so that conformational variability is best seen
%             implemented only for entities with a single chain
%             if there is no folded domain, superposition is in inertia
%             frame
%
%   entity = ADAPTIVE_SUPERPOSITION(entity,options,eau)
%
% INPUT
% entity    entity that describes the ensemble, must be single peptide
%           chain
% options   computation options, struct with fields
%           .threshold  ensemble aligned uncertainty threshold for folded
%                       domains, defaukts to 12.5 A
%           .minsize    minimum size of a folded domain, defaults to 25
%           .minlink    minimum length of a flexible linker between
%                       domains, defaults to 3
%           .interact   threshold for relative ensemble aligned uncertainty
%                       for domain interaction, defaults to 5/3
%           .ensemble   optional name for ensemble file list of
%                       superimposed ensemble, if missing or empty, no
%                       ensemble file is written
%           .minconf    minimum number of conformers in a subensemble,
%                       defaults to 3
% eau       ensemble aligned uncertainty matrix, optional, is computed when
%           missing, can also be provided as field of entity, if the
%           argument is present, it overwrites eau information in entity
%
% OUTPUT
% entity       entity after superposition, coordinates have changed, domain
%              definition is added or updated, ensemble aligned uncertainty
%              is added whend missing

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2025: Gunnar Jeschke

options.eau = true;
options.rbfile = '';
options.script = '';
options.unify = 15;
options.minterm = 10;
options.local = 5;
options.figname ='';

if ~isfield(options,'threshold') || isempty(options.threshold)
    options.threshold = 12.5;
end

if ~isfield(options,'minsize') || isempty(options.minsize)
    options.minsize = 25;
end

if ~isfield(options,'minlink') || isempty(options.minlink)
    options.minlink = 3;
end

if ~isfield(options,'interact') || isempty(options.interact)
    options.interact = 5/3;
end

if ~isfield(options,'minconf')  || isempty(options.minconf)
    options.minconf = 3;
end

if exist('eau','var')
    entity.eau = eau;
end

if ~isfield(entity,'eau')
    entity = ensemble_aligned_uncertainty(entity);
end

[R,~] = size(entity.eau);
entity.eau_resnum = 1:R;
entity.eau_chain = 'A';

entity = domain_partitioning(entity,options);

entity = inertia_frame(entity);

[ndom,~] = size(entity.domains);
deau = zeros(1,ndom);
sortings = cell(1,ndom);
all_rmsd = zeros(1,ndom);
all_diff = zeros(1,ndom);
C = length(entity.populations);
if ndom > 0
    merits = zeros(1,ndom);
    assignments = zeros(ndom,C);
    for dom = 1:ndom
        ceau = entity.eau(entity.domains(dom,1):entity.domains(dom,2),entity.domains(dom,1):entity.domains(dom,2));
        deau(dom) = mean(ceau,'all')';
        fprintf(1,'Mean EAU for domain %i is %5.2f Å\n',dom,deau(dom));
        selected = sprintf('(%s)%i-%i',entity.eau_chain,...
            entity.eau_resnum(entity.domains(dom,1)),...
            entity.eau_resnum(entity.domains(dom,2)));
        [centity,rmsd] = superimpose_ensemble(entity,selected);
        all_rmsd(dom) = sqrt(sum(rmsd.^2/length(rmsd)));
        fprintf(1,'Superposition rmsd for domain %i is %5.2f Å\n',dom,all_rmsd(dom));
        D = pair_rmsd_matrix_oriented(centity,entity.eau_chain);
        [C,~] = size(D);
        weights = zeros(C);
        for c1 = 1:C-1
            for c2 = c1+1:C
                weights(c1,c2) = abs(c1-c2);
                weights(c2,c1) = weights(c1,c2);
            end
        end
        weights = weights.^2;
        diff = 1/sum(sum(weights.*D));
        indices = 1:C;
        maxpasses = 150;
        passes = 0;
        poi = 0;
        changed = true;
        while passes < maxpasses && changed
            passes = passes +1;
            changed = false;
            for c1 = 1:C-1
                for c2 = c1+1:C
                    poi = poi + 1;
                    indices2 = indices;
                    indices2(c1) = indices(c2);
                    indices2(c2) = indices(c1);
                    D2 = D(indices2,indices2);
                    dtest = 1/sum(sum(weights.*D2));
                    if dtest < diff
                        diff = dtest;
                        indices = indices2;
                        changed = true;
                    end
                end
            end
        end        
        fprintf(1,'Sorting quality for superposition upon domain %i is %12.6g\n',dom,diff);
        all_diff(dom) = diff;
        sortings{dom} = indices; 
        figure; hold on;
        image(D(indices,indices),'CDataMapping','scaled');
        [assignment,merit] = ensemble_partionining(D(indices,indices),options);
        merits(dom) = merit;
        assignments(dom,:) = assignment;
        fprintf(1,'Merit of partitioning into %i ensembles: %6.4f\n',max(assignment),merit);
        curr_axis = gca;
        curr_axis.YDir = 'normal';
        c = colorbar;
        c.Label.String = 'RMSD [Å]';
        axis tight
        xlabel('Conformer');
        ylabel('Conformer');
        axis equal
        title(sprintf('Superposition on domain %i',dom));
    end
    if max(merits) <= 1
        [mindiff,bdom] = min(all_diff);
        fprintf(1,'Best domain for superposition is %i at %12.6g Å\n',bdom,mindiff);
    else
        [merit,bdom] = max(merits);
        assignment = assignments(bdom,:);
        fprintf(1,'Best domain for superposition is %i at subsensemble merit %6.4f with %i subensembles\n',bdom,merit,max(assignment));
    end
    selected = sprintf('(%s)%i-%i',entity.eau_chain,...
        entity.eau_resnum(entity.domains(bdom,1)),...
        entity.eau_resnum(entity.domains(bdom,2)));
    centity = superimpose_ensemble(entity,selected);
    indices = sortings{bdom};
else % no folded domain detected
    centity = entity;
    D = pair_drms_matrix(centity,entity.eau_chain);
    [C,~] = size(D);
    weights = zeros(C);
    for c1 = 1:C-1
        for c2 = c1+1:C
            weights(c1,c2) = abs(c1-c2);
            weights(c2,c1) = weights(c1,c2);
        end
    end
    weights = weights.^2;
    diff = 1/sum(sum(weights.*D));
    indices = 1:C;
    maxpasses = 150;
    passes = 0;
    all_diff = zeros(1,maxpasses*C*(C-1)/2);
    poi = 0;
    changed = true;
    while passes < maxpasses && changed
        passes = passes +1;
        changed = false;
        for c1 = 1:C-1
            for c2 = c1+1:C
                poi = poi + 1;
                indices2 = indices;
                indices2(c1) = indices(c2);
                indices2(c2) = indices(c1);
                D2 = D(indices2,indices2);
                dtest = 1/sum(sum(weights.*D2));
                all_diff(poi) = dtest;
                if dtest < diff
                    diff = dtest;
                    indices = indices2;
                    changed = true;
                end
            end
        end
    end
    figure; hold on;
    image(D(indices,indices),'CDataMapping','scaled');
    [assignment,merit] = ensemble_partionining(D(indices,indices),options);
    curr_axis = gca;
    curr_axis.YDir = 'normal';
    c = colorbar;
    c.Label.String = 'RMSD [Å]';
    axis tight
    xlabel('Conformer');
    ylabel('Conformer');
    axis equal
    title('dRMSD sorting');
    fprintf(1,'Found %i subensembles with merit %6.4f\n',max(assignment),merit);
end
entity.subensembles = assignment;

if isfield(options,'ensemble') && ~isempty(options.ensemble)
    if ~contains(options.ensemble,'.')
        options.ensemble = strcat(options.ensemble,'.ens');
    end
    ens_fid = fopen(options.ensemble,'wt');
    fprintf(ens_fid,'%% Sorted and superimposed ensemble %s\n',centity.name);
    for conf = 1:length(centity.populations)
        fname = sprintf('%s_m%i.pdb',centity.name,conf);
        clear save_options
        save_options.order = indices(conf);
        put_pdb(centity,fname,save_options);
        fprintf(ens_fid,'%s  %8.6f\n',fname,centity.populations(indices(conf)));
    end
end
fclose(ens_fid);

function [indices,merit] = ensemble_partionining(D,options)

if ~exist('options','var') || ~isfield(options,'minconf')
    options.minconf = 3;
end

Z = linkage(D,'complete');

[n,~] = size(D);
all_indices = zeros(5,n);
all_merits = zeros(1,5);

for s = 2:5
    indices = cluster(Z,'maxclust',s);
    own = 0;
    other = 0;
    valid = true;
    for k = 1:s
        se_size = sum(indices==k);
        if se_size < options.minconf
            valid = false;
        end
        own = own + sqrt(mean(D(indices==k,indices==k).^2,"all"));
    end
    own = own/s;
    for k1 = 1:s-1
        for k2 = k1+1:s
            other = other + sqrt(mean(D(indices==k1,indices==k2).^2,"all"));
        end
    end
    other = 2*other/(s*(s-1));
    if valid
        merit = other/own;
    else
        merit = 0;
    end
    all_indices(s,:) = indices;
    all_merits(s) = merit;
end
[merit,best] = max(all_merits);
if merit < 1
    merit = 0;
    indices = ones(1,n);
else
    indices = all_indices(best,:);
end