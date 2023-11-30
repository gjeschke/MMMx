function [results,PAE_distribution] = AF2_human_proteome_analysis(params,interactive)

cap = 31.75; % cap of AlphaFold2 PAE

if ~exist('interactive','var') || isempty(interactive)
    interactive = false;
end

unify = 15;

DB{50000} = ''; % allocate space for cell array
nprot = 0;
ifid = fopen('AF2_human_proteome.txt');
sfid = fopen('AF2_proteome_sequences.fasta','wt');
while 1
    tline = fgetl(ifid);
    if ~ischar(tline) 
        break 
    end
	nprot = nprot + 1;
	DB{nprot} = tline;
end
fclose(ifid);
DB = DB(1:nprot);

fid = fopen('AF2_human_proteome_analysis.log','wt');

results.F = zeros(length(params.minsize),1); % completely folded
results.Sf = zeros(length(params.minsize),1); % single domain/flexible terminus
results.Mf = zeros(length(params.minsize),1); % multiple domains/flexible terminus
results.D = zeros(length(params.minsize),1); % fully disordered

lengths = zeros(1,nprot);

n_valid = 0;

long_PAE_axis = 0:0.5:32;
mean_long_PAE = zeros(1,length(long_PAE_axis));

for nr = 1:nprot

    if mod(nr,10) == 0
        fprintf(1,'%4.1f completed\n',100*(nr-1)/nprot);
    end

    options = weboptions('ContentType','json');
    basname = DB{nr};
    poi = strfind(basname,'-model');
    parts = split(basname,'-');

    fname = sprintf('%s-predicted_aligned_error_v4.json',basname(1:poi-1));
    url = sprintf('https://alphafold.ebi.ac.uk/files/%s',fname);

    try
        val = webread(url,options);
    catch
        fprintf(fid,'For protein %i the predicted alignment error file could not be read.\n\n',nr);
        continue
    end

    pairdist = (val.predicted_aligned_error + val.predicted_aligned_error')/2;
    pairdist0 = pairdist;

    [n,~] = size(pairdist);

    if n < 50
        fprintf(fid,'Protein %i has a sequence length of only %i residues and is skipped.\n\n',n);
        continue
    end

    UniProtID = parts{2};
    query = sprintf('https://alphafold.ebi.ac.uk/api/prediction/%s',UniProtID);
    try
        AF_info = webread(query);
    catch
        continue;
    end
    fprintf(sfid,'>sp|%s\n',AF_info.uniprotAccession);
    tsequence = AF_info.uniprotSequence;
    while length(tsequence)>70
        fprintf(sfid,'%s\n',tsequence(1:70));
        tsequence=tsequence(71:end);
    end
    fprintf(sfid,'%s\n\n',tsequence);
    if interactive
        figure(1); clf; hold on
        image(pairdist,'CDataMapping','scaled');
        curr_axis = gca;
        set(curr_axis,'YDir','normal');
        colorbar;
        axis tight
        xlabel('Residue number');
        ylabel('Residue number');
        axis equal
        title(sprintf('Protein:%i (%s)',nr,DB{nr}));
    end

    [n,~] = size(pairdist);
    
    num = 0;
    mean_PAE = 0;
    for k1 = 1:n-10
        for k2 = k1+10:n
            num = num + 1;
            mean_PAE = mean_PAE + pairdist(k1,k2);
        end
    end
    mean_PAE = mean_PAE/num;
    [~,index] = min(abs(mean_PAE-long_PAE_axis));
    mean_long_PAE(index) = mean_long_PAE(index) + 1; 

    n_valid = n_valid + 1;
    lengths(nr) = n;

    for parset = 1:length(params.minsize)

        minsize = params.minsize(parset);
        mindomain = params.mindomain(parset);
        minlink = params.minlink(parset);
        minterminal = params.minterminal(parset);
        threshold = params.threshold(parset);
        fprintf(fid,'Analysis of the AlphaFold2 human proteome with the following parameters:');
        fprintf(fid,'Minimum protein size: %i\n',minsize);
        fprintf(fid,'Minimum domain size: %i\n',mindomain);
        fprintf(fid,'Minimum interdomain linker size: %i\n',minlink);
        fprintf(fid,'Minimum length of flexible terminus: %i\n',minterminal);
        fprintf(fid,'Folded domain threshold: %5.1f%c\n\n',threshold,char(197));
        
        pairdist = pairdist0;

        A = pairdist < threshold;

        local = 5;

        extension = ones(1,n);
        for k1 = 2:n-1
            ext = 1;
            k2 = k1+local;
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
        b = ones(1,local)/local;
        a = 1;
        extension = filter(b,a,extension);
        
        dpoi = 0;
        domains = zeros(ceil(n/mindomain),2);
        while max(extension) >= mindomain
            [ext,k] = max(extension);
            dstart = k;
            dend = k + round(ext);
            if dend > n
                dend = n;
            end
            mean_uncert0 = mean(mean(pairdist(dstart:dend,dstart:dend)));
            for k1 = dstart:-1:1
                mean_uncert = mean(pairdist(k1,k1:dend));
                if mean_uncert <= threshold
                    dstart = k1;
                    mean_uncert0 = mean(mean(pairdist(dstart:dend,dstart:dend)));
                else
                    break;
                end
            end
            if mean_uncert0 > threshold
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
                if mean_uncert <= threshold
                    mean_uncert0 = mean_uncert;
                    dend = k1;
                else
                    break;
                end
            end
            if mean_uncert0 > threshold
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
            for k = 1:n
                test = k + round(extension(k));
                if test <= n && extension(test) > cap
                    extension(k) = 1;
                end
            end
            if mean_uncert0 <= threshold
                dpoi = dpoi +1;
                domains(dpoi,:) = [dstart,dend];
            end
        end
        
        % unify domains separated by too short linkers, if requested
        if minlink > 1
            cdpoi = 1;
            for k = 2:dpoi
                if domains(k,1) - domains(k-1,2) < minlink + 1
                    test = max(max(pairdist(domains(k-1,1):domains(k,2),domains(k-1,1):domains(k,2))));
                    if test <= unify
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
        if max(max(domains)) == 0
            dpoi = 0;
        else
            % sort domains
            [~,indices] = sort(domains(1:dpoi,1));
            domains = domains(indices,:);
        end

        
        folded = 0;
        for k = 1:dpoi
            k1 = domains(k,1);
            k2 = domains(k,2);
            folded = folded + k2 - k1 + 1;
            if interactive
                plot([k1,k1],[k1,k2],'LineWidth',2,'Color',[0.8,0,0]);
                plot([k2,k2],[k1,k2],'LineWidth',2,'Color',[0.8,0,0]);
                plot([k1,k2],[k1,k1],'LineWidth',2,'Color',[0.8,0,0]);
                plot([k1,k2],[k2,k2],'LineWidth',2,'Color',[0.8,0,0]);
            end
        end
%         if interactive
%             keyboard
%         end
        
        max_flex_link = 0;
        overlap = 0;
        for k = 2:dpoi
            if domains(k,1)-domains(k-1,2) - 1 > max_flex_link
                max_flex_link = domains(k,1)-domains(k-1,2);
            end
            if domains(k,1) < domains(k-1,2)
                overlap = 1;
            end
        end
        
        max_terminal = NaN;
        if dpoi > 0
            max_terminal = domains(1,1) - 1;
            if n - domains(dpoi,2) > max_terminal
                max_terminal = n - domains(end,2);
            end
        end
        
        % classify this protein
        if dpoi < 1
            results.D(parset) = results.D(parset) + 1;
        else
            if dpoi > 1
                results.Mf(parset) = results.Mf(parset) + 1;
            elseif max_terminal == 0
                results.F(parset) = results.F(parset) + 1;
            elseif max_terminal < minterminal
                results.F(parset) = results.F(parset) + 1;
            else
                results.Sf(parset) = results.Sf(parset) + 1;
            end
        end
        
        if overlap
            fprintf(2,'Overlapping domains detected for protein %i\n',nr);
        end
        if interactive
            drawnow
        end
       
    end
    
end

fclose(fid);
fclose(sfid);

for parset = 1:length(params.minsize)
    normalizer = sum(results.F(parset)+results.Sf(parset)+results.Mf(parset)+results.D(parset));
    results.F(parset) = results.F(parset)/normalizer;
    results.Sf(parset) = results.Sf(parset)/normalizer;
    results.Mf(parset) = results.Mf(parset)/normalizer;
    results.D(parset) = results.D(parset)/normalizer;
end

PAE_distribution.axis = long_PAE_axis;
PAE_distribution.count = mean_long_PAE;
