clear entities
clear options
options.fname = 'FUS_sorted2_acs.csv';
options.visualization = 'FUS_sorted2.mmm';
addresses = {};
entities{1} = get_PED('PED00494','e001');

[coor,coor_3D,measures,D] = abstract_conformer_space(entities,addresses,options);

save FUS_sorted2_ensemble_acs coor coor_3D measures D

figure(111);
C = length(measures.convergence.rmean);
plot(1:C,measures.convergence.rmean,'.','MarkerSize',12);

figure(112);
plot(1:C,measures.convergence.stdr,'.','MarkerSize',12);

figure(113);
plot(1:C,measures.convergence.Rg_acs,'.','MarkerSize',12);

figure(114);
plot(1:C,measures.convergence.Rg_real,'.','MarkerSize',12);

figure(115);
plot(1:C,measures.convergence.dimension,'.','MarkerSize',12);

figure(116);
plot(1:C,measures.convergence.resolution,'.','MarkerSize',12);

L_curve_x = (1:C)/C;
L_curve_y = measures.convergence.dimension/max(measures.convergence.dimension);
[~,natdim2] = min(L_curve_x.^2+L_curve_y.^2);
fprintf(1,'The natural dimension by the corner method is %i with embedding error of %4.1f %c\n',...
    natdim2,tinv(0.975,C-1)*measures.convergence.dimension(natdim2),char(197));

fprintf(1,'The 3D embedding error is %4.1f %c\n',measures.error_3D,char(197));
fprintf(1,'ACS resolution is %4.1f %c\n',measures.resolution,char(197));
fprintf(1,'3D coordinates by %s\n',measures.mode_3D);
fprintf(1,'Disorder ratio %5.3f\n',measures.Rg_acs/measures.Rg);
