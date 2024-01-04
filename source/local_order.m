function local_order(entity,fname,figures)
%
% LOCAL_ORDER Computes an order parameter for a residue with respect to the
%               complete structure
%
%   LOCAL_ORDER(entity,fname)
%   Save results in files based on fname, if empty, only plots are
%   displayed
%
%   LOCAL_ORDER(entity,fname,figures)
%   Save plots in given figure format
%
%   based on a reinterpretation of Flory's characteristic ratio for
%   polymers that may feature defined secondary structure
%
% INPUT
% entity    entity in MMMx:atomic format 
% fname     basis file name for output, extension is disregared and
%           replaced with .csv, data is save as comma-separated value file,
%           one file for each chain <fname>_<cid>.csv, where <cid> is the
%           chain identifier
% figures   specifies figure format for saving as anb extension, if empty 
%           or missing, figures are not saved, options are:
%           fig (Matlab figure), jpg, png, eps, pdf, bmp, emf, pbm, pcx,
%           pgm, ppm, svg, tif
%           recommended: 'pdf' or 'svg' for use with professional software
%                        'emf' for use with Microsoft Word
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% maximum number of atoms for array pre-allocation, function gets slow, if
% this number is too small and is memory-intensive, if it is too large

disordered.mean = 0.1906;
disordered.std = 0.0007;


% initialize empty outputs
if ~exist('fname','var')
    fname = '';
end

if ~isempty(fname) % remove file extension if any
    [fpath,fnamex,~] = fileparts(fname);
    fname = fullfile(fpath,fnamex);
end

if ~exist('figures','var')
    figures = '';
elseif ~isempty(figures) && figures(1) == '.' % be tolerant if user specified extension with leading dot
    figures = figures(2:end);
end

if isempty(fname)
    figname0 = sprintf('local_order_%s',entity.name);
else
    figname0 = fname;
end

chains = fieldnames(entity);
pop = entity.populations/sum(entity.populations); % you can never normalize too often
for c = 1:length(chains)
    if isstrprop(chains{c},'upper') % this is a chain field
        chain = chains{c};
        residues = fieldnames(entity.(chain));
        all_monomers = cell(1,length(pop));
        for co = 1:length(pop)
            monomers = zeros(length(residues),3);
            resnums = zeros(length(residues),1);
            respoi = 0;
            for r = 1:length(residues)
                residue = residues{r};
                if residue(1) == 'R' % this is a residue field
                    resnum = str2double(residue(2:end));
                    type = 'hetero';
                    if isfield(entity.(chain).(residue),'N') % if an N atom record exists, we assume an amino acid
                        type = 'aa';
                    end
                    if isfield(entity.(chain).(residue),'P') % if an P atom record exists, we assume an nucleic acid
                        type = 'na';
                    end
                    switch type
                        case 'aa'
                            if ~isfield(entity.(chain).(residue),'CA') % if CA atom of this residue is missing, phi is undetermined
                                continue
                            end
                            respoi = respoi + 1;
                            resnums(respoi) = resnum;
                            index = entity.(chain).(residue).CA.tab_indices(co);
                            monomers(respoi,:) = entity.xyz(index,:);
                        case 'na'
                            if ~isfield(entity.(chain).(residue),'C4_') % if C4' atom of this nucleic acid is missing, phi is undetermined
                                continue
                            end
                            respoi = respoi + 1;
                            resnums(respoi) = resnum;
                            index = entity.(chain).(residue).C4_.tab_indices(co);
                            monomers(respoi,:) = entity.xyz(index,:);
                    end
                end
            end
            monomers = monomers(1:respoi,:);
            all_monomers{co} = monomers;
        end
        resnums = resnums(1:respoi);
        order_poi = 0;
        order = zeros(respoi-1,1);
        res_order = zeros(respoi-1,1);
        for i = 1:respoi-1
            if resnums(i+1) ~= resnums(i) + 1 % the "bond" vector for this residue must be defined
                continue
            end
            order_poi = order_poi + 1;
            res_order(order_poi) = resnums(i);
            corr_cosij = zeros(1,respoi-1);
            ij_poi = 0;
            for j = 1:respoi-1
                if resnums(j+1) ~= resnums(j) + 1 % the "bond" vector for the second residue must be defined
                    continue
                end
                ij_poi = ij_poi + 1;
                all_cosij = zeros(size(pop));
                for co = 1:length(pop)
                    monomers = all_monomers{co};
                    bond_i = monomers(i+1,:) - monomers(i,:);
                    bond_i = bond_i/norm(bond_i);
                    bond_j = monomers(j+1,:) - monomers(j,:);
                    bond_j = bond_j/norm(bond_j);
                    all_cosij(co) = sum(bond_i.*bond_j); % cosine is the scalar product of the normalized vectors
                end
                mean_cosij = sum(pop.*all_cosij);
                var_cosij = sum(pop.*(all_cosij - mean_cosij).^2);
                corr_cosij(ij_poi) = 1 - sqrt(2*var_cosij);
            end
            corr_cosij = corr_cosij(1:ij_poi);
            order(order_poi) = sum(corr_cosij)/ij_poi;
        end
        order = order(1:order_poi);
        res_order = res_order(1:order_poi);
        h = figure; clf; hold on;
        ub = disordered.mean + 2*disordered.std;
        lb = disordered.mean - 2*disordered.std;
        liner = ones(size(res_order'));
        fill([res_order', fliplr(res_order')],[ub*liner, lb*liner],0.75*[1,1,1],'LineStyle','none');
        plot(res_order,disordered.mean*liner,'-','LineWidth',2,'Color',0.25*[1,1,1]);
        plot(res_order,order,'.','MarkerSize',14,'Color',[0,0,0.7]);
        xlabel('Residue number');
        ylabel('Site-specific order');
        mi = min(order) - 0.05;
        if mi > disordered.mean - 2*disordered.std-0.05
            mi = disordered.mean - 2*disordered.std-0.05;
        end
        if mi < 0
            mi = 0;
        end
        ma = max(order) + 0.05;
        if ma < disordered.mean + 2*disordered.std + 0.05
            ma = disordered.mean + 2*disordered.std + 0.05;
        end
        if ma > 1
            ma = 1;
        end
        axis([0,max(res_order)+1,mi,ma]);
        title(sprintf('%s chain %s',entity.name, chain));
        if ~isempty(figures)
            figname = sprintf('%s_chain_%s.%s',figname0,chain,figures);
            saveas(h,figname);
        end
        if ~isempty(fname)
            outmat = [res_order order];
            datname = sprintf('%s_chain_%s.csv',fname,chain);
            writematrix(outmat,datname)
        end
    end
end