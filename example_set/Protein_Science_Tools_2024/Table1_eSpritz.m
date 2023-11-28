mindomains = 25;
thresholds = 0.14:0.001:0.15;

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

results = AF2_human_proteome_analysis_eSpritz(params);
normalize = mean(results.F+results.Sf+results.Mf+results.D);

save AF2_human_proteome_analysis_results params results
fprintf(1,'Folded: (%4.1f +/- %4.1f%%)\n',100*mean(results.F)/normalize,200*std(results.F));
fprintf(1,'Single domain, flexible termini: (%4.1f +/- %4.1f%%)\n',100*mean(results.Sf)/normalize,200*std(results.Sf));
fprintf(1,'Multi-domain: (%4.1f +/- %4.1f%%)\n',100*mean(results.Mf)/normalize,200*std(results.Mf));
fprintf(1,'Fully disordered: (%4.1f +/- %4.1f%%)\n',100*mean(results.D)/normalize,200*std(results.D));

