clear
tlines{1} = '23, A, 1, 92'; % native, NMR, SAXS
tlines{2} = '159, A, 1, 92'; % native, NMR, SAXS
tlines{3} = '160, A, 1, 92'; % native, NMR, SAXS, smFRET
tlines{4} = '424, A, 1, 92'; % native, IDPConformerGenerator, chem. shift data
tlines{5} = '423, A, 1, 92'; % native, IDPConformerGenerator, no exp. data
tlines{6} = '455, A, 1, 92'; % native, idpGAN machine learning
tlines{7} = '487, A, 1, 92'; % native, idpGAN machine learning
tlines{8} = '1, A, 1, 92'; % phosphorylated 5, 33, 45, 69, 76, 80, NMR, SAXS
tlines{9} = '486, A, 1, 92'; % phosphorylated 5, 33, 45, 69, 76, 80, idpGAN machine learning
tlines{10} = '161, A, 1, 92'; % phosphorylated 2, 5, 33, 45, 69, 76, 80, NMR, SAXS, smFRET
tlines{11} = '454, A, 1, 92'; % phosphorylated 2, 5, 33, 45, 69, 76, 80, idpGAN machine learning
tlines{12} = '14, A, 1, 92'; % phosphorylated 5, 33, 45, 69, 76, 80, complex with SPK1 and CDC4

ens = 0;
ensembles(100).entity = [];
ensembles(100).chain = '';
ensembles(100).range = [];
for l = 1:length(tlines)
    args = split(tlines{l},',');
    pedID = strtrim(args{1});
    while length(pedID) < 5
        pedID = strcat('0',pedID);
    end
    pedID = strcat('PED',pedID);
    chain = strtrim(args{2});
    if length(args) >= 4
        fr = str2double(args{3});
        lr = str2double(args{4});
        range = [round(fr),round(lr)];
    else
        range = [];
    end
    fprintf(1,'Importing %s for analyzing (%s)%i-%i\n',pedID,chain,range);
    entities = get_PED(pedID);
    if ~iscell(entities)
        entity = entities;
        clear entities
        entities{1} = entity;
    end
    for k = 1:length(entities)
        ens = ens + 1;
        ensembles(ens).entity = entities{k};
        ensembles(ens).chain = chain;
        ensembles(ens).range = range;
        ensembles(ens).pedID = pedID;
    end
end
ensembles = ensembles(1:ens);

save Sic1_full_set ensembles

similarities = ones(ens);
Rg = ones(ens,1);
poi = 0;
for ens1 = 1:length(ensembles)-1
    chain1 = ensembles(ens1).chain;
    range1 = ensembles(ens1).range;
    entity1 = ensembles(ens1).entity;
    for ens2 = ens1+1:length(ensembles)
        poi = poi + 1;
        fprintf(1,'Comparing %i,%i (%5.1f%%)\n',ens1,ens2,100*poi/(ens*(ens-1)/2));
        chain2 = ensembles(ens2).chain;
        range2 = ensembles(ens2).range;
        entity2 = ensembles(ens2).entity;
        [similarity,~,~,~,Rg1,Rg2] = ensemble_similarity(entity1,chain1,range1,entity2,chain2,range2);
        similarities(ens1,ens2) = similarity;
        similarities(ens2,ens1) = similarity;
        Rg(ens1) = Rg1;
        Rg(ens2) = Rg2;
    end
end

save Sic1_full_set ensembles similarities Rg

figure(1); clf; hold on
my_map = flipud(colormap);
colormap(my_map);
image(similarities,'CDataMapping','scaled');
curr_axis = gca;
set(curr_axis,'YDir','normal');
colorbar;
axis tight
xlabel('Ensemble number');
ylabel('Ensemble number');
axis equal
title('Comparison of Sic1 ensembles');
set(gca,'FontSize',12);
