cap = 31.75;

mindomains = 25;
thresholds = cap/3-2:0.1:cap/3+2;

ntests = length(mindomains)*length(thresholds);

params.minsize = 50*ones(ntests,1);
params.minlink = 3*ones(ntests,1);
params.minterminal = 10*ones(ntests,1);
params.mindomain = zeros(ntests,1);
params.threshold = zeros(ntests,1);

poi = 0;
for k1 = 1:length(mindomains)
    for k2 = 1:length(thresholds)
        poi = poi + 1;
        params.mindomain(poi) = mindomains(k1); 
        params.threshold(poi) = thresholds(k2); 
    end
end

[results,PAE_distribution] = AF2_human_proteome_analysis(params,false);

save AF2_human_proteome_analysis_results params results PAE_distribution

figure(1); clf; plot(thresholds,results.F);
figure(2); clf; plot(thresholds,results.Sf);
figure(3); clf; plot(thresholds,results.Mf);
figure(4); clf; plot(thresholds,results.D);

fprintf(1,'Folded: (%4.1f +/- %4.1f%%)\n',100*mean(results.F),200*std(results.F));
fprintf(1,'Single domain, flexible termini: (%4.1f +/- %4.1f%%)\n',100*mean(results.Sf),200*std(results.Sf));
fprintf(1,'Multi-domain: (%4.1f +/- %4.1f%%)\n',100*mean(results.Mf),200*std(results.Mf));
fprintf(1,'Fully disordered: (%4.1f +/- %4.1f%%)\n',100*mean(results.D),200*std(results.D));

figure(5);
bar(PAE_distribution.axis,PAE_distribution.count);
axis([min(PAE_distribution.axis),max(PAE_distribution.axis),0,1.05*max(PAE_distribution.count)]);
set(gca,'FontSize',12);
xlabel('PAE');
ylabel('Protein count');
