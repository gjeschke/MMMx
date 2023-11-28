function results = AF2_human_proteome_analysis_eSpritz(params,interactive)

if ~exist('interactive','var') || isempty(interactive)
    interactive = false;
end

DB{50000} = ''; % allocate space for cell array
nprot = 0;
ifid = fopen('AF2_human_proteome.txt');
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

fid = fopen('AF2_human_proteome_analysis_eSpritz.log','wt');

results.F = zeros(length(params.minsize),1); % completely folded
results.Sf = zeros(length(params.minsize),1); % single domain/flexible terminus
results.Mf = zeros(length(params.minsize),1); % multiple domains/flexible terminus
results.D = zeros(length(params.minsize),1); % fully disordered

n_valid = 0;

for nr = 1:nprot

    if mod(nr,10) == 0
        fprintf(1,'%4.1f completed\n',100*(nr-1)/nprot);
    end
    
    n_valid = n_valid + 1;

    basname = DB{nr};
    parts = split(basname,'-');
    UniProtID = parts{2};


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

    [~,probability] = get_disorder_espritz(UniProtID);
    if isempty(probability)
        continue
    end

    for parset = 1:length(params.minsize)

        threshold = params.threshold(parset);
        mindomain = params.mindomain(parset);
        minlink = params.minlink(parset);
        minterminal = params.minterminal(parset);
        disorder = probability > threshold;
        [~,domains,max_terminal] = classify_by_disorder(disorder,minterminal,mindomain,minlink);
        fprintf(fid,'Analysis of the AlphaFold2 human proteome with the following parameters:');
        fprintf(fid,'Minimum domain size: %i\n',mindomain);
        fprintf(fid,'Minimum interdomain linker size: %i\n',minlink);
        fprintf(fid,'Minimum length of flexible terminus: %i\n',minterminal);
        fprintf(fid,'Folded domain threshold: %5.1f%c\n\n',threshold,char(197));

        [dpoi,~] = size(domains);
        
        D = 0;
        M = 0;
        F = 0;
        Fs = 0;
        Ft = 0;
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

results.F = results.F/n_valid;
results.Sf = results.Sf/n_valid;
results.Mf = results.Mf/n_valid;
results.D = results.D/n_valid;

