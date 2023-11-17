function [classifier,domains,max_terminal] = classify_by_disorder(disorder,minterminal,mindomain,minlink)

if ~exist('minterminal','var')
    minterminal = 10;
end

if ~exist('mindomain','var')
    mindomain = 25;
end

if ~exist('minlink','var')
    minlink = 3;
end


% maximum length of flexible N-terminal section
poi = 0;
while length(disorder) <= poi + 1 && disorder(poi+1) == 1
    poi = poi + 1;
end
max_terminal = poi;
res = length(disorder) + 1;

% maximum length of flexible C-terminal section
poi = 0;
while res - poi - 1 > 0 && disorder(res-poi-1) == 1
    poi = poi + 1;
end
if poi > max_terminal
    max_terminal = poi;
end

% find ordered domains
dpoi = 0;
domains = zeros(1000,2);
dstart = 0;
dend0 = 0;
while dstart < length(disorder)
    while dstart < length(disorder) && disorder(dstart+1) == 1 
        dstart = dstart + 1;
    end
    if dstart < length(disorder)
        dend = dstart;
        while dend < length(disorder) && disorder(dend+1) == 0
            dend = dend + 1;
        end
    else
        continue
    end
    if dend-dstart+1 >= mindomain
        dpoi = dpoi + 1;
        domains(dpoi,:) = [dstart,dend];
    end
    if dend ~= dend0
        dstart = dend+1;
        dend0 = dend;
    end
end

% unify domains separated by too short linkers, if requested
if minlink > 1
    cdpoi = 1;
    for k = 2:dpoi
        if domains(k,1) - domains(k-1,2) < minlink + 1
                domains(cdpoi,2) = domains(k,2);
        else
            cdpoi = cdpoi + 1;
        end
    end
    dpoi = cdpoi;
end

if max(max(domains)) == 0
    dpoi = 0;
end

domains = domains(1:dpoi,:);

if dpoi < 1
    classifier = 'D';
else
    if dpoi > 1
        classifier = 'M';
    elseif max_terminal == 0
        classifier = 'F';
    elseif max_terminal < minterminal
        classifier = 'Fs';
    else
        classifier = 'Ft';
    end
end
