function pdb_parcels(entity,chain_id,basname)
% PDB_PARCELS(entity,basname)
%
% Distributes an entity over several PDB files so that each PDB file has
% less than 100'000 atoms, with a maximum of 26 conformers per file,  
% conformers are saved as coordinate-separated individual chains
% only one chain of the conformer is selected
%
% Input:
%
% entity    entity specifying the whole ensemble
% chain_id  chain identifier of the chain to be analyzed
% basname   basis file name, individual files are <basname>_p%i.pdb, whhere
%           %i is the file number
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2023: Gunnar Jeschke

% Transform coordinates to dind a compact representation
entity = inertia_frame(entity,sprintf('(%c)',chain_id));

nres = 0;
residues = fieldnames(entity.(chain_id));
for kr = length(residues) % expand over all residues
    residue = residues{kr};
    if strcmp(residue(1),'R') % these are residue fields
        resnum = str2double(residue(2:end));
        if resnum > nres
            nres = resnum;
        end
    end
end

% determine box sides
ax = 10 + max(entity.xyz(:,1)) - min(entity.xyz(:,1));
ay = 10 + max(entity.xyz(:,2)) - min(entity.xyz(:,2));
az = 10 + max(entity.xyz(:,3)) - min(entity.xyz(:,3));

C = length(entity.populations); % number of conformers
na = length(entity.elements)/C; % number of atoms per conformer

cf = floor(99999/na); % number of conformers per file
if cf > floor(9999/nres)
    cf = floor(9999/nres);
end
if cf > 27
    cf = 27;
end

nf = ceil(C/cf); % number of files
c = 0; % initialize conformer number

for kf = 1:nf % loop over individual PDB files
    fname = sprintf('%s_p%i.pdb',basname,kf);
    kc = 0; % conformer index in this PDB file
    options.fid = -1;
    options.datnum = 0;
    options.dresnum = 0;
    while kc < cf && c < C
        kc = kc + 1;
        c = c + 1;
        kz = floor(kc/9);
        ky = floor((kc-9*kz)/3);
        kx = kc-9*kz-3*ky;
%         chain = char(double('A')-1+kc);
%         options.chainIDs{1,1} = chain_id;
%         options.chainIDs{1,2} = chain;
        % compute coordinate offset
        dx = ax * (kx-1);
        dy = ay * (ky-1);
        dz = az * (kz-1);
        options.dxyz = [dx,dy,dz]; % set coordinate offset
        options.selected = true;
        entity.selected = c;
        if kc < cf && c < C
            [~,fid] = put_pdb(entity,fname,options);
            options.fid = fid;
        else
            put_pdb(entity,fname,options);
        end
        options.datnum = options.datnum + na;
        options.dresnum = options.dresnum + nres;
    end
end
