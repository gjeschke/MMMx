function outdata = flexibility_analysis(entity,fname,figures)
%
% FLEXIBILITY_ANALYSIS Analyze Ramachandran angle distributions in
%                      ensemble, for RNA, pseudo-torsion angles are analyzed 
%
%   FLEXIBILITY_ANALYSIS(entity,fname)
%   Save results in files based on fname, if empty, only plots are
%   displayed
%
%   FLEXIBILITY_ANALYSIS(entity,fname,figures)
%   Save plots in given figure format
%
%   based on circular variance, see: 
%   http://www.fiserlab.org/manuals/procheck/manual/man_cv.html
%
% INPUT
% entity    entity in MMMx:atomic format 
% fname     basis file name for output, extension is disregared and
%           replaced with .csv, data is saved as comma-separated value file,
%           one file for each chain <fname>_<cid>.csv, where <cid> is the
%           chain identifier
% figures   specifies figure format for saving as an extension, if empty 
%           or missing, figures are not saved, options are:
%           fig (Matlab figure), jpg, png, eps, pdf, bmp, emf, pbm, pcx,
%           pgm, ppm, svg, tif
%           recommended: 'pdf' or 'svg' for use with professional software
%                        'emf' for use with Microsoft Word
%
% OUPUT
% outdata   cell(1,chains) of matrices (residues,3) with residue numbers in
%           column 1, site flexibility in column 2, and residue type in
%           column 3, residue types are: 1 for non-special residues, 2 for
%           Gly, 3 for Pro, and 4 for Thr

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

% maximum number of atoms for array pre-allocation, function gets slow, if
% this number is too small and is memory-intensive, if it is too large

% Flexibility parameters for residue types from analysis of unrestrained 
% FUS NTD
flexpar(2).mean = 0.7158; % GLY
flexpar(2).std = 0.0226;
flexpar(3).mean = 0.2574; % PRO
flexpar(3).std = 0.0089;
flexpar(4).mean = 0.3437; % THR
flexpar(4).std = 0.0051;
flexpar(1).mean = 0.4065; % all others
flexpar(1).std = 0.0203;


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
    figname0 = sprintf('flexibility_%s',entity.name);
else
    figname0 = fname;
end

chains = fieldnames(entity);
outdata = cell(1,length(chains));
chain_count = 0;
pop = entity.populations/sum(entity.populations); % you can never normalize too often
for c = 1:length(chains)
    if isstrprop(chains{c},'upper') % this is a chain field
        chain = chains{c};
        chain_count = chain_count + 1;
        residues = fieldnames(entity.(chain));
        flexibilities = zeros(length(residues),1);
        resnums = zeros(length(residues),1);
        restypes = ones(length(residues),1);
        respoi = 0;
        for r = 1:length(residues)
            residue = residues{r};
            if residue(1) == 'R' % this is a residue field
                resnum = str2double(residue(2:end));
                previous = sprintf('R%i',resnum-1);
                next = sprintf('R%i',resnum+1);
                if ~isfield(entity.(chain),previous) % previous residue must exist
                    continue
                end
                if ~isfield(entity.(chain),next) % next residue must exist
                    continue
                end
                type = 'hetero';
                if isfield(entity.(chain).(residue),'N') % if an N atom record exists, we assume an amino acid
                    type = 'aa';
                end
                if isfield(entity.(chain).(residue),'P') % if an P atom record exists, we assume an nucleic acid
                    type = 'na';
                end
                switch type
                    case 'aa'
                        if ~isfield(entity.(chain).(previous),'C') % if C atom of previous residue is missing, phi is undetermined
                            continue
                        end
                        if ~isfield(entity.(chain).(residue),'CA') % if CA atom of this residue is missing, phi is undetermined
                            continue
                        end
                        if ~isfield(entity.(chain).(residue),'C') % if C atom of this residue is missing, phi is undetermined
                            continue
                        end
                        if ~isfield(entity.(chain).(next),'N') % if N atom of next residue is missing, psi is undetermined
                            continue
                        end
                        phi = zeros(length(pop),1);
                        psi = zeros(length(pop),1);
                        for co = 1:length(pop)
                            index = entity.(chain).(previous).C.tab_indices(co);
                            Cim = entity.xyz(index,:); % coordinate of previous C atom
                            index = entity.(chain).(residue).N.tab_indices(co);
                            Ni = entity.xyz(index,:); % coordinate of N atom
                            index = entity.(chain).(residue).CA.tab_indices(co);
                            CAi = entity.xyz(index,:); % coordinate of CA atom
                            index = entity.(chain).(residue).C.tab_indices(co);
                            Ci = entity.xyz(index,:); % coordinate of C atom
                            phi(co) = dihedral_mmmx(Cim,Ni,CAi,Ci);
                            index = entity.(chain).(next).N.tab_indices(co);
                            Nip = entity.xyz(index,:); % coordinate of N atom
                            psi(co) = dihedral_mmmx(Ni,CAi,Ci,Nip);                                                        
                        end
                        Rav = sqrt((sum(cos(phi))^2 + sum(sin(phi))^2 + sum(cos(psi))^2 + sum(sin(psi))^2)/2)/length(phi);
                        respoi = respoi + 1;
                        resnums(respoi) = resnum;
                        flexibilities(respoi) = 1-Rav;
                        restype = upper(entity.(chain).(residue).name);
                        switch restype
                            case 'GLY'
                                restypes(respoi) = 2;
                            case 'PRO'
                                restypes(respoi) = 3;
                            case 'THR'
                                restypes(respoi) = 4;
                        end
                    case 'na'
                        if ~isfield(entity.(chain).(previous),'C4_') % if C4' atom of previous nucleic acid is missing, phi is undetermined
                            continue
                        end
                        if ~isfield(entity.(chain).(residue),'C4_') % if C4' atom of this nucleic acid is missing, phi is undetermined
                            continue
                        end
                        if ~isfield(entity.(chain).(next),'P') % if P atom of next nucleic acid is missing, psi is undetermined
                            continue
                        end
                        eta = zeros(length(pop),1); % needed for proper circular mean
                        theta = zeros(length(pop),1);
                        for co = 1:length(pop)
                            index = entity.(chain).(previous).C4_.tab_indices(co);
                            C4pm = entity.xyz(index,:); % coordinate of previous C atom
                            index = entity.(chain).(residue).P.tab_indices(co);
                            P = entity.xyz(index,:); % coordinate of N atom
                            index = entity.(chain).(residue).C4_.tab_indices(co);
                            C4p = entity.xyz(index,:); % coordinate of CA atom
                            index = entity.(chain).(next).P.tab_indices(co);
                            Pp = entity.xyz(index,:); % coordinate of P atom next na
                            index = entity.(chain).(next).C4_.tab_indices(co);
                            C4pp = entity.xyz(index,:); % coordinate of C4 atom'next na
                            eta(co) = dihedral_mmmx(C4pm,P,C4p,Pp);
                            theta(co) = dihedral_mmmx(P,C4p,Pp,C4pp);                                                        
                        end
                        Rav = sqrt((sum(pop.*cos(eta))^2 + sum(pop.*sin(eta))^2 + sum(pop.*cos(theta))^2 + sum(pop.*sin(theta))^2)/2);
                        respoi = respoi + 1;
                        resnums(respoi) = resnum;
                        flexibilities(respoi) = 1-Rav;
                end
            end
        end
        resnums = resnums(1:respoi);
        flexibilities = flexibilities(1:respoi);
        restypes = restypes(1:respoi);
        h = figure; clf; hold on;
        colors = [0.25,0.25,0.25;0.7,0,0;0,0.5,0;0,0,0.7];
        for res = 1:length(resnums)
            errorbar(resnums(res),flexpar(restypes(res)).mean,2*flexpar(restypes(res)).std,...
                'Color',colors(restypes(res),:)+[0.3,0.3,0.3]);
            plot(resnums(res),flexibilities(res),'.','MarkerSize',14,'Color',colors(restypes(res),:));
        end
        xlabel('Residue number');
        ylabel('Site-specific flexibility');
        title(sprintf('%s chain %s',entity.name, chain));
        axis([0,resnums(end)+1,0,1]);
        if ~isempty(figures)
            figname = sprintf('%s_chain_%s.%s',figname0,chain,figures);
            saveas(h,figname);
        end
        outmat = [resnums flexibilities restypes];
        outdata{chain_count} = outmat;
        if ~isempty(fname)
            datname = sprintf('%s_chain_%s.csv',fname,chain);
            writematrix(outmat,datname)
        end
    end
end

outdata = outdata(1:chain_count);