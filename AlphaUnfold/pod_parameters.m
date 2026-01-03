function [f_IDR,f_fuzzy,f_residual,nd] = pod_parameters(UniProtID)

options.structure = false;
options.pLDDT = false;
entity = get_AF(UniProtID,options);
psize = length(entity.AF_info.uniprotSequence);

f_IDR = entity.AF_info.fractionPlddtVeryLow + entity.AF_info.fractionPlddtLow;
if floor(f_IDR*psize) > 15 % at least 15 disordered residues
    f_residual = entity.AF_info.fractionPlddtLow/...
        (entity.AF_info.fractionPlddtVeryLow + entity.AF_info.fractionPlddtLow);
else
    f_residual = 1;
end

if floor((1-f_IDR)*psize) > 15 % at least 15 residues in IFRs
    f_fuzzy = entity.AF_info.fractionPlddtConfident/(entity.AF_info.fractionPlddtVeryHigh + entity.AF_info.fractionPlddtConfident);
else
    f_fuzzy = 0;
end

domains = get_domains(entity.pae);
[nd,~] = size(domains);
