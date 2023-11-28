function [results,PAE_distribution] = AF2_human_proteome_analysis(params,interactive)

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

        % classes of proteins
        F = 0;
        Ft = 0;
        Fs = 0;
        M = 0;
        D = 0;
        
        A = pairdist < threshold;

        local = 5;

        for k1 = 1:n-2*local
            for k2 = k1-local:k1+local
                if k2 >=1 && k2 <=n
                    A(k1,k2) = 0;
                    A(k2,k1) = 0;
                end
            end
        end

        % find diagonal blocks, based on F. Pedroche, M. Rebollo, C. Carrascosa and A. Palomares (2012)
        % http://arxiv.org/abs/1206.5726
        L = diag( sum(A,2) ) - A;
        value = sum( triu(L) );
        freq = find( value == 0 );
        
        dpoi = 0;
        domains = zeros(length(freq),2);
        dstart = 1;
        if length(freq) ==1
            mean_uncert = mean(mean(pairdist));
            if mean_uncert <= threshold
                dpoi = dpoi + 1;
                domains(dpoi,1) = dstart;
                domains(dpoi,2) = freq(1);
            end
        else
            for k = 2:length(freq)
                if freq(k)-dstart >= minsize
                    mean_uncert = mean(mean(pairdist(dstart:freq(k),dstart:freq(k))));
                    if mean_uncert <= threshold
                        dpoi = dpoi + 1;
                        domains(dpoi,1) = dstart;
                        domains(dpoi,2) = freq(k);
                    end
                end
                dstart = freq(k);
            end
        end
        
        % refine domain boundaries by sliding window approach
        
        for dom = 1:dpoi
            dstart = domains(dom,1);
            dend = domains(dom,2);
            ka = dstart - 10;
            if ka < 1
                ka = 1;
            end
            ke = dstart + 10;
            if ke > n
                ke = n;
            end
            determ = zeros(1,ke-ka+1);
            for k = ka:ke
                determ(k-ka+1) = (mean(pairdist(k,dstart:dend))+mean(pairdist(dstart:dend,k)))/2;
            end
            [~,k] = min(abs(determ - threshold));
            domains(dom,1) = ka+k-1;
            ka = dend - 10;
            if ka < 1
                ka = 1;
            end
            ke = dend + 10;
            if ke > n
                ke = n;
            end
            determ = zeros(1,ke-ka+1);
            for k = ka:ke
                determ(k-ka+1) = (mean(pairdist(k,dstart:dend))+mean(pairdist(dstart:dend,k)))/2;
            end
            [~,k] = min(abs(determ - threshold));
            domains(dom,2) = ka+k-1;
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
            D = D + 1;
        else
            if dpoi > 1
                M = M + 1;
            elseif max_terminal == 0
                F = F + 1;
            elseif max_terminal < minterminal
                Fs = Fs + 1;
            else
                Ft = Ft + 1;
            end
        end
        
        if overlap
            fprintf(2,'Overlapping domains detected for protein %i\n',nr);
        end
        if interactive
            drawnow
        end
       
        results.F(parset) = results.F(parset) + F + Fs;
        results.Sf(parset) = results.Sf(parset) + Ft;
        results.Mf(parset) = results.Mf(parset) + M;
        results.D(parset) = results.D(parset) + D;

    end
    
end

fclose(fid);
fclose(sfid);

results.F = results.F/n_valid;
results.Sf = results.Sf/n_valid;
results.Mf = results.Mf/n_valid;
results.D = results.D/n_valid;

PAE_distribution.axis = long_PAE_axis;
PAE_distribution.count = mean_long_PAE;
