function [paevec,seqdist,audata,IFR_pairs,type] = get_pae_from_UniProt(UniProtID,options)

paevec = [];
seqdist = [];
audata = [];
IFR_pairs = [];
type = [];

cap = 31; % cap of AlphaFold2 PAE

if ~exist('options','var') || ~isfield(options,'threshold') || isempty(options.threshold)
    options.threshold = cap/3;
end

if ~isfield(options,'unfify') || isempty(options.unify)
    options.unify = 15;
end

if ~isfield(options,'minsize') || isempty(options.minsize)
    options.minsize = 25;
end

if ~isfield(options,'minlink') || isempty(options.minlink)
    options.minlink = 3;
end

if ~isfield(options,'local') || isempty(options.local)
    options.local = 5;
end

if ~isfield(options,'minterm') || isempty(options.minterm)
    options.minterminal = 10;
end

if ~isfield(options,'termini') || isempty(options.termini)
    options.termini = true;
end

if ~isfield(options,'tight_threshold') || isempty(options.tight_threshold)
    options.tight_threshold = 5;
end

if ~isfield(options,'interactive') || isempty(options.interactive)
    options.interactive = false;
end

entity = get_AF(UniProtID,options);

if isempty(entity) || ~isfield(entity,'pae') || isempty(entity.pae)
    return
end

pae = entity.pae;
paevec = pae(:);

[n,~] = size(pae);

seqdist = repmat(1:n, n, 1) - repmat((1:n)', 1, n);
seqdist = seqdist(:);

pairdist = (pae + pae')/2;

A = pairdist < options.threshold;

extension = ones(1,n);
for k1 = 2:n-1
    ext = 1;
    k2 = k1+options.local;
    while k2 < n
        k2 = k2 + 1;
        if A(k1,k2)
            ext = k2-k1;
        else
            break
        end
    end
    extension(k1) = ext;
end
% moving average filter
b = ones(1,options.local)/options.local;
a = 1;
extension = filter(b,a,extension);

ndom = 0;
domains = zeros(ceil(n/options.minsize),2);
while max(extension) >= options.minsize
    [ext,k] = max(extension);
    dstart = k;
    dend = k + round(ext);
    if dend > n
        dend = n;
    end
    mean_uncert0 = mean(mean(pairdist(dstart:dend,dstart:dend)));
    for k1 = dstart:-1:1
        mean_uncert = mean(pairdist(k1,k1:dend));
        if mean_uncert <= options.threshold
            dstart = k1;
            mean_uncert0 = mean(mean(pairdist(dstart:dend,dstart:dend)));
        else
            break;
        end
    end
    if mean_uncert0 > options.threshold
        for k1 = dstart:n
            mean_uncert = mean(mean(pairdist(k1:dend,k1:dend)));
            if mean_uncert < mean_uncert0
                mean_uncert0 = mean_uncert;
                dstart = k1;
            else
                break;
            end
        end
    end
    for k1 = dend:n
        mean_uncert = mean(mean(pairdist(k1,dstart:k1)));
        if mean_uncert <=options.threshold
            mean_uncert0 = mean_uncert;
            dend = k1;
        else
            break;
        end
    end
    if mean_uncert0 > options.threshold
        for k1 = dend:-1:1
            mean_uncert = mean(mean(pairdist(dstart:k1,dstart:k1)));
            if mean_uncert < mean_uncert0
                mean_uncert0 = mean_uncert;
                dend = k1;
            else
                break;
            end
        end
    end
    extension(dstart:dend) = zeros(1,dend-dstart+1);
    dompae = pairdist(dstart:dend,dstart:dend);
    pairdist(dstart:dend,dstart:dend) = 1e12;
    for k = 1:n
        test = k + round(extension(k));
        if test <= n && extension(test) > cap
            extension(k) = 1;
        end
    end
    if mean_uncert0 <= options.tight_threshold && max(max(dompae)) < cap
        ndom = ndom +1;
        domains(ndom,:) = [dstart,dend];
    end
end

domains = domains(1:ndom,:);

switch ndom
    case 0
        type = 4;
    case 1
        if domains(1,1) > options.minterminal || n-domains(1,2) > options.minterminal
            type = 2;
        else
            type = 1;
        end
    otherwise
        type = 3;
end

audata = zeros(ndom-1,2);
IFR_pairs = 0;
dpae = cell(1,ndom);
for dom = 1:ndom
    dpae{dom} = pae(domains(dom,1):domains(dom,2),domains(dom,1):domains(dom,2));
end
for dom = 1:ndom-1
    ipae1 = pae(domains(dom,1):domains(dom,2),domains(dom+1,1):domains(dom+1,2));
    ipae2 = pae(domains(dom+1,1):domains(dom+1,2),domains(dom,1):domains(dom,2));
    if max(max(ipae1)) >= cap || max(max(ipae2)) >= cap
        continue
    else
        IFR_pairs = IFR_pairs + 1;
        audata(IFR_pairs,1) = domains(dom,2) - domains(dom,1);
        ipae = mean((ipae1+ipae2')/2,'all');
        audata(IFR_pairs,2) = round(ipae);
    end
end
