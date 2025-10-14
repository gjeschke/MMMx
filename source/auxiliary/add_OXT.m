function entity = add_OXT(entity,chain)
% entity = add_OXT(entity,chain)
%
% Adds a terminal oxygen to a chain made by the Flex algorithm
%
% Input:
%
% entity    MMMx entity
% chain     char, chain identifier, usually 'A'
%
% Output:
%
% entity    modified MMMx entity
%
% G. Jeschke, 9.10.2025

b = 1.261; % C-OXT bond length
bvec = [b,0,0]; % C-OXT bond vector along x
fi = 120*pi/180; % C_alpha-C-OXT angle
sfi = sin(fi); 
cfi = cos(fi);


% Determine the last residue
residues = fieldnames(entity.(chain));
maxresnum = 0;
for kr = 1:length(residues)
    residue = residues{kr};
    if strcmp(residue(1),'R') % these are residue fields
        resnum = str2double(residue(2:end));
        if resnum > maxresnum
            maxresnum = resnum;
            restag = residue;
        end
    end
end

if maxresnum > 0 % entity is changed only if the chain has at least one residue
    if ~isfield(entity.(chain).(restag),'OXT')
        % local frame of Calpha of the C-terminal residue
        CA_index = entity.(chain).(restag).CA.tab_indices(1);
        N_index = entity.(chain).(restag).N.tab_indices(1);
        C_index = entity.(chain).(restag).C.tab_indices(1);
        CA = entity.xyz(CA_index,:);
        N = entity.xyz(N_index,:);
        C = entity.xyz(C_index,:);
        x = C - CA;
        x = x/norm(x);
        yp = CA - N;
        yp = yp/norm(yp);
        z = cross_rowvec(x,yp);
        z = z/norm(z);
        y = cross_rowvec(z,x);
        y = y/norm(y);
        A = [x;y;z];
        % transformation matrix into that local frame
        A = A';
        ctau = cos(pi);
        stau = sin(pi);
        A2 = [-cfi,-sfi,0;sfi*ctau,-cfi*ctau,-stau;sfi*stau,-cfi*stau,ctau];
        A = A*A2;
        OXT = C + A*bvec'; % coordinates of OXT
        xyz = entity.xyz;
        [n,~] = size(xyz);
        xyz = [xyz;OXT];
        entity.xyz = xyz;
        entity.(chain).(restag).OXT.element = 'O';
        entity.(chain).(restag).OXT.charge = 0;
        entity.(chain).(restag).OXT.bfactor = 0;
        entity.(chain).(restag).OXT.selected = 0;
        entity.(chain).(restag).OXT.selected_locations = 1;
        entity.(chain).(restag).OXT.tab_indices = n+1;
        % determine last atom index
        atoms = fieldnames(entity.(chain).(restag));
        atindex = 0;
        for ka = 1:length(atoms)
            atom = atoms{ka};
            if isstrprop(atom(1),'upper') % these are atom fields
                atindex = atindex + 1;
            end
        end
        entity.(chain).(restag).OXT.index = atindex + 1;
        entity.occupancies = [entity.occupancies;100];
    end
end