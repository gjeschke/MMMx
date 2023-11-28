mindomains = 25;
thresholds = 10;

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

[results,PAE_distribution] = AF2_human_proteome_analysis(params);

figure(1);
bar(PAE_distribution.axis,PAE_distribution.count);
axis([min(PAE_distribution.axis),max(PAE_distribution.axis),0,1.05*max(PAE_distribution.count)]);
set(gca,'FontSize',12);
xlabel('PAE');
ylabel('Protein count');
