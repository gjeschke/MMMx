datafile = 'A1_minus_buffer_concnorm_cut.dat';
pdbfile = 'shnRNPA1_conformers_32_231_c1';

options.lm = 20;
options.fb = 18;
options.delete = true;
options.cst = false;
options.err = true;
options.crysol3 = true;

[fit,chi2,outname,status,result] = fit_SAXS(datafile,pdbfile,options);

fprintf(1,'chi^2 = %6.3f\n',chi2);

figure(1); clf;
plot(fit(:,1),fit(:,2),'k');
hold on;
plot(fit(:,1),fit(:,3),'r');
if options.err
    plot(fit(:,1),fit(:,4),'b');
end


