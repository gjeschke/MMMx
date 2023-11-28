multi_ensemble = 0;
multi_ensembles(500).ID = ''; 

tic,
for id = 1:500
    pedID = sprintf('%i',id);
    while length(pedID) < 5
        pedID = strcat('0',pedID);
    end
    query = sprintf('https://deposition.proteinensemble.org/api/v1/entries/PED%s',pedID);
    try
        PED_info = webread(query);
    catch exception
        exceptions{1} = exception;
        continue
    end
    last_id = id;
    if length(PED_info.ensembles) > 1
        feasible = true;
        for ens = 1:length(PED_info.ensembles)
            if PED_info.ensembles(1).models > 1000
                feasible = false;
            end
        end
        if feasible
            multi_ensemble = multi_ensemble + 1;
            multi_ensembles(multi_ensemble).ID = pedID;
            multi_ensembles(multi_ensemble).ensembles = PED_info.ensembles;
        end
    end
end
toc,

fprintf(1,'%i PED entries were checked\n',last_id);
fprintf(1,'%i PED entries have multiple ensembles\n',multi_ensemble);
multi_ensembles = multi_ensembles(1:multi_ensemble);

tic,
for k = 1:multi_ensemble
    n = length(multi_ensembles(k).ensembles);
    ensembles = cell(1,n);
    for ens = 1:n
        entity = get_PED(sprintf('PED%s',multi_ensembles(k).ID),multi_ensembles(k).ensembles(ens).ensemble_id);
        ensembles{ens} = entity;
    end
    fprintf(1,'%i ensembles read for PED entry %s. Now comparing...\n',n,multi_ensembles(k).ID);
    multi_ensembles(k).similarities = zeros(n);
    multi_ensembles(k).Rgdet = zeros(n);
    for kn1 = 1:n-1
        for kn2 = kn1+1:n
            ens2 = ensembles{kn2};
            [similarity,~,~,~,Rg1,Rg2] = ensemble_similarity(ensembles{kn1},'',[],ensembles{kn2},'',[]);
            Rg_dist = sqrt((Rg1-Rg2)^2/(Rg1*Rg2));
            multi_ensembles(k).Rgdet(kn1,kn2) = Rg_dist;
            multi_ensembles(k).Rgdet(kn2,kn1) = Rg_dist;
            multi_ensembles(k).similarities(kn1,kn2) = similarity;
            multi_ensembles(k).similarities(kn2,kn1) = similarity;
        end
    end
end
toc,

save PED_multiensemble_similarities multi_ensembles

data = load('PED_multiensemble_similarities');

n = length(data.multi_ensembles);

distribution = zeros(1,101);
sim_axis = linspace(0,1,length(distribution));

all_similarities = zeros(1,1000);
all_Rgdet = zeros(1,1000);
all_sizes = zeros(1,1000);
others = distribution;
tauf = distribution;
mvn = distribution;
pairs = 0;
for k = 1:n
    similarities = data.multi_ensembles(k).similarities;
    [ens,~] = size(similarities);
    min_sim = 1;
    for e1 = 1:ens-1
        for e2 = e1+1:ens
            pairs = pairs + 1;
            all_similarities(pairs) = similarities(e1,e2);
            ID = str2double(data.multi_ensembles(k).ID);
            if all_similarities(pairs) < 0.9
                fprintf(1,'Entry %s, ensembles (%i,%i) are dissimilar with %6.4f\n',data.multi_ensembles(k).ID,e1,e2,all_similarities(pairs));
            end
            all_Rgdet(pairs) = data.multi_ensembles(k).Rgdet(e1,e2);
            all_sizes(pairs) = data.multi_ensembles(k).ensembles(e1).models*data.multi_ensembles(k).ensembles(e2).models;
            poi = 1 + round(100*similarities(e1,e2));
            distribution(poi) = distribution(poi) + 1;
            switch ID
                case 17
                    tauf(poi) = tauf(poi) + 1;
                case 20
                    mvn(poi) = mvn(poi) + 1;
                otherwise
                    others(poi) = others(poi) + 1;
            end
            if similarities(e1,e2) < min_sim
                min_sim = similarities(e1,e2);
            end
        end
    end
end

all_similarities = all_similarities(1:pairs);
all_Rgdet = all_Rgdet(1:pairs);
all_sizes = all_sizes(1:pairs);

pairs = sum(distribution);
fprintf(1,'\n%i ensemble pairs were compared\n',pairs);

figure(1); clf;
bar(sim_axis,distribution);
axis([0,1,0,1.05*max(distribution)]);
set(gca,'FontSize',12);
xlabel('Similarity s_{kl}');
ylabel('Number of ensemble pairs');

figure(2); clf;
b = bar(sim_axis,[others;tauf;mvn],'stacked');
h = [b(2),b(3),b(1)];
axis([0,1,0,1.05*max(distribution)]);
set(gca,'FontSize',12);
xlabel('Similarity s_{kl}');
ylabel('Number of ensemble pairs');
legend(h,'PED00017','PED00020','Others','Location','northwest')