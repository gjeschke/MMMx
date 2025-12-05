function domains = get_domains(pae,options)

cap = 31;

if ~exist('options','var') ||~isfield(options,'threshold') || isempty(options.threshold)
    options.threshold = cap/3;
end

if ~isfield(options,'unify') || isempty(options.unify)
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
 if ~isfield(options,'interact') || isempty(options.interact)
    options.interact = 1;
end
pairdist = (pae + pae')/2;
[n,~] = size(pairdist);

A = pairdist < options.threshold;
pairdist0 = pairdist;
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

dpoi = 0;
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
    % for k1 = dend:n
    %     mean_uncert = mean(mean(pairdist(k1,dstart:k1)));
    %     if mean_uncert <=options.threshold
    %         mean_uncert0 = mean_uncert;
    %         dend = k1;
    %     else
    %         break;
    %     end
    % end
    span = dstart - 1;
    if n - dend > span
        span = n - dend;
    end
    dstart0 = dstart;
    dend0 = dend;
    for k = 1:span
        extended = false;
        if dend+k <= n
            mean_uncert = mean(mean(pairdist(dend0+k,dstart:dend0+k)));
            if mean_uncert <= options.threshold
                mean_uncert0 = mean_uncert;
                dend = dend0 + k;
                extended = true;
            end
        end
        if dstart-k >= 1
            mean_uncert = mean(mean(pairdist(dstart0-k,dstart0-k:dend)));
            if mean_uncert <= options.threshold
                mean_uncert0 = mean_uncert;
                dstart = dstart0 - k;
                extended = true;
            end
        end
        if ~extended
            break
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
    pairdist(dstart:dend,dstart:dend) = 1e12;
    % for k = 1:n
    %     test = k + round(extension(k));
    %     if test <= n && extension(test) > cap
    %         extension(k) = 1;
    %     end
    % end
    if mean_uncert0 <= options.threshold && dend - dstart >= options.minsize 
        dpoi = dpoi +1;
        domains(dpoi,:) = [dstart,dend];
    end
end

% sort domains
[~,indices] = sort(domains(1:dpoi,1));
domains = domains(indices,:);

pairdist = pairdist0;
% unify domains separated by too short linkers, if requested
if options.minlink > 1
    cdpoi = 1;
    for k = 2:dpoi
        if domains(k,1) - domains(k-1,2) < options.minlink + 1
            test = max(max(pairdist(domains(k-1,1):domains(k,2),domains(k-1,1):domains(k,2))));
            if test <= options.unify
                domains(cdpoi,2) = domains(k,2);
            elseif domains(k,1) - domains(k-1,2) < 3
                domains(k,1) = domains(k,1) + 1;
                domains(k-1,2) = domains(k-1,2) - 1;
                cdpoi = cdpoi + 1;
            else
                cdpoi = cdpoi + 1;
            end
        else
            cdpoi = cdpoi + 1;
        end
    end
    dpoi = cdpoi;
end

% unify domains that are strongly interacting

combine = true;
while combine
    combine = false;
    cdpoi = 1;
    for k = 2:dpoi
        mpae1 = mean(pairdist(domains(k-1,1):domains(k-1,2),domains(k-1,1):domains(k-1,2)),"all");
        mpae2 = mean(pairdist(domains(k,1):domains(k,2),domains(k,1):domains(k,2)),"all");
        mpaec = mean(pairdist(domains(k-1,1):domains(k,2),domains(k-1,1):domains(k,2)),"all");
        defn = 2*mpaec/(mpae1+mpae2);
        if defn <= options.interact
            combine = true;
            domains(cdpoi,2) = domains(k,2);
        else
            cdpoi = cdpoi + 1;
        end
    end
    dpoi = cdpoi;
end

if isempty(domains) || max(max(domains)) == 0
    dpoi = 0;
else
    % sort domains
    [~,indices] = sort(domains(1:dpoi,1));
    domains = domains(indices,:);
end

domains = domains(1:dpoi,:);
