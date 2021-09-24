function f = ensemble_distance(pair_rmsd,popA,popB)

distinct = 1e-3; % threshold for considering ensembles as distinct

% cardinalities
CA = length(popA);
CB = length(popB);


SAB = 0;
for cA = 1:CA
    for cB = 1:CB
        SAB = SAB + popA(cA)*popB(cB)*pair_rmsd(cA,CA+cB);
    end
end

SA = 0;
for c1 = 1:CA
    for c2 = 1:CA
        SA = SA + popA(c1)*popA(c2)*pair_rmsd(c1,c2);
    end
end

SB = 0;
for c1 = 1:CB
    for c2 = 1:CB
        SB = SB + popB(c1)*popB(c2)*pair_rmsd(c1+CA,c2+CA);
    end
end

f = SAB - (SA + SB)/2;

if f < distinct
    f = 0;
end