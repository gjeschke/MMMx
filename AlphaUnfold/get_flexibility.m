function flexibility = get_flexibility(entity,GPR_disorder,GPR_RMSF,GPR_SSP)
% flexibility = get_flexibility(entity,GPR_disorder,GPR_RMSF,GPR_SSP)
%     Extracts protein order-disorder parameters for a protein contained in
%     the Big Fantastic Virus Database

flexibility.fIDR = [];
flexibility.ffuzzy = [];
flexibility.fresidual = [];
flexibility.pLDDT = [];
flexibility.nd = [];
flexibility.disorder = [];
flexibility.RMSF = [];
flexibility.SSP = [];

pLDDT = entity.pLDDT;
psize = length(pLDDT);
flexibility.pLDDT = pLDDT;


fractionPlddtVeryLow = sum(pLDDT < 50)/psize;
fractionPlddtLow = sum((pLDDT >= 50) & (pLDDT < 70))/psize;
fractionPlddtConfident = sum((pLDDT >= 70) & (pLDDT < 90))/psize;
fractionPlddtVeryHigh = sum(pLDDT >= 90)/psize;

fIDR = fractionPlddtVeryLow + fractionPlddtLow;
flexibility.fIDR = fIDR;
if floor(fIDR*psize) > 15 % at least 15 disordered residues
    flexibility.fresidual = fractionPlddtLow/(fractionPlddtVeryLow + fractionPlddtLow);
else
    flexibility.fresidual = 1;
end
if floor((1-fIDR)*psize) > 15 % at least 15 residues in IFRs
    flexibility.ffuzzy = fractionPlddtConfident/(fractionPlddtVeryHigh + fractionPlddtConfident);
else
    flexibility.ffuzzy = 0;
end

cap = 31;

domain_options.threshold = 0.65*cap;
domain_options.unify = 15;
domain_options.minsize = 15;
domain_options.minlink = 3;
domain_options.local = 15;
domain_options.interact = 1;

domains = get_domains(entity.pae,domain_options);
[flexibility.nd,~] = size(domains);

disorder = get_disorder(pLDDT,entity.pae,GPR_disorder);
flexibility.disorder = sum(disorder)/length(disorder);
RMSF = get_RMSF(pLDDT,entity.pae,GPR_RMSF);
flexibility.RMSF = mean(RMSF(disorder == 0));
SSP = get_SSP(pLDDT,entity.pae,GPR_SSP);
flexibility.SSP = mean(SSP(disorder ~= 0));

