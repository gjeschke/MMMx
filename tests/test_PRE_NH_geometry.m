entity = get_pdb('1nti');

NHbond = 0.98;

rNH = zeros(1,85);
ang_err = zeros(1,85);
H_dev = zeros(1,85);
poi = 0;
for res = 2:86
    resfield = sprintf('R%i',res);
    resfield_prev = sprintf('R%i',res-1);
    if isfield(entity.A.(resfield),'H')
        poi = poi + 1;
        H_index = entity.A.(resfield).H.tab_indices(1);
        N_index = entity.A.(resfield).N.tab_indices(1);
        CA_index = entity.A.(resfield).CA.tab_indices(1);
        C_index = entity.A.(resfield_prev).C.tab_indices(1);
        HN = entity.xyz(H_index,:)-entity.xyz(N_index,:);
        rNH(poi) = norm(HN);
        CN = entity.xyz(C_index,:)-entity.xyz(N_index,:);
        CAN = entity.xyz(CA_index,:)-entity.xyz(N_index,:);
        dir_HN_exp = -(CN/norm(CN)+CAN/norm(CAN));
        dir_HN_exp = dir_HN_exp/norm(dir_HN_exp);
        dir_HN = HN/norm(HN);
        H_pred = entity.xyz(N_index,:) + NHbond*dir_HN_exp;
        H_true = entity.xyz(H_index,:);
        H_dev(poi) = norm(H_pred - H_true);
        ang_err(poi) = 180*acos(sum(dir_HN_exp.*dir_HN))/pi;
    end
end


rNH = rNH(1:poi);
ang_err = ang_err(1:poi);
H_dev = H_dev(1:poi);

figure(1); clf;
plot(rNH,'k.');
title('H-N bond length');

figure(2); clf;
plot(ang_err,'k.');
title('Angle deviation');

figure(3); clf;
plot(H_dev,'k.');
title('Deviation of proton position');


