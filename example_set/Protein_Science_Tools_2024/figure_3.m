panel = 'a';
switch panel
    case 'a'
        UniProtID = 'A6PVI3';
    case 'b'
        UniProtID = 'Q6IPF2';
    case 'c'
        UniProtID = 'P26599';
    case 'd'
        UniProtID = 'A6NN06';
end
rigiflex_options.label = 'mtsl';
rigiflex_options.threshold = 10;

entity = get_AF(UniProtID);
entity = domain_partitioning(entity,rigiflex_options);

[disorder_eSpritz,probability] = get_disorder_espritz(UniProtID);
[disorder_SETH,CheZOD] = get_disorder_SETH('classification_test_set_SETH_disorder.fasta',UniProtID);


[classifier,eSpritz_domains] = classify_by_disorder(disorder_eSpritz);
fprintf(1,'eSpritz classifies %s as %s\n',UniProtID,classifier);


[classifier,SETH_domains] = classify_by_disorder(disorder_SETH);
fprintf(1,'SETH classifies %s as %s\n',UniProtID,classifier);

[dpoi,~] = size(eSpritz_domains); 
for k = 1:dpoi
    k1 = eSpritz_domains(k,1);
    k2 = eSpritz_domains(k,2);
    plot([k1,k1],[k1,k2],'LineWidth',2,'Color',[1,0.48,0.63]);
    plot([k2,k2],[k1,k2],'LineWidth',2,'Color',[1,0.48,0.63]);
    plot([k1,k2],[k1,k1],'LineWidth',2,'Color',[1,0.48,0.63]);
    plot([k1,k2],[k2,k2],'LineWidth',2,'Color',[1,0.48,0.63]);
end

[dpoi,~] = size(SETH_domains); 
for k = 1:dpoi
    k1 = SETH_domains(k,1);
    k2 = SETH_domains(k,2);
    plot([k1,k1],[k1,k2],'LineWidth',2,'Color',[1,1,1]);
    plot([k2,k2],[k1,k2],'LineWidth',2,'Color',[1,1,1]);
    plot([k1,k2],[k1,k1],'LineWidth',2,'Color',[1,1,1]);
    plot([k1,k2],[k2,k2],'LineWidth',2,'Color',[1,1,1]);
end

set(gca,'FontSize',12);

figure; hold on;
plot(disorder_eSpritz);
[dpoi,~] = size(eSpritz_domains);
for k = 1:dpoi
    plot(eSpritz_domains(k,:),[0.05,0.05],'LineWidth',2,'Color',[0.75,0,0]);
end
figure; hold on;
plot(disorder_SETH);
[dpoi,~] = size(SETH_domains);
for k = 1:dpoi
    plot(SETH_domains(k,:),[0.05,0.05],'LineWidth',2,'Color',[0.75,0,0]);
end
