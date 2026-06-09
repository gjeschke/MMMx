function [disorder,scores] = get_disorder(pLDDT,pae,GPR)

X = zeros(length(pLDDT),length(GPR.feature_selection));

for r = 1:length(pLDDT)
    resi = r - GPR.frame_width;
    if resi < 1
        resi = 1;
    end
    rese = r + GPR.frame_width;
    if rese > length(pLDDT)
        rese = length(pLDDT);
    end
    this_pae = pae(resi:rese,resi:rese);
    this_pLDDT = pLDDT(resi:rese);
    features = extract_features(this_pae,this_pLDDT);
    features = features(GPR.feature_selection);
    X(r,:) = features;
end
scores = predict(GPR.best_model, X);
disorder = scores > GPR.level;
