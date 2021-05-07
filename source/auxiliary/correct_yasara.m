function exceptions = correct_yasara(fname,oname,substitute,newcid)
% function exception = correct_yasara(fname,oname,substitute,newcid)
% 
% Corrects non-standard atom names in Yasara files and optionally
% substitutes chain identifiers
%
% Input:
%
% fname         input PDB file
% oname         output PDB file, defaults to [fname '_corr.pdb']
% substitute    string of chain identifiers that are to be substituted
% newcid        string of new chain identifiers for substitution
%
% Output:
%
% exceptions    cell(1,1) of a potential Matlab exception, if input or
%               output file could not be opened
%
% G. Jeschke, 7.5.2021

exceptions{1} = [];

if ~exist('substitute','var')
    substitute = '';
end
if ~exist('newcid','var')
    newcid = substitute;
end

if ~contains(fname,'.pdb')
    fname = strcat(fname,'.pdb');
end
ifid=fopen(fname);
if ifid==-1
    exceptions = {MException('correct_yasara:input_file_could_not_be_opened',...
                                    'Input file %s could not be opened',fname)};
    return;
end

if exist('oname','var') && ~isempty(oname)
    ofid=fopen(oname,'wt');
else
    ofid=fopen(strcat(fname,'_corr.pdb'),'wt');
end
if ofid==-1
    exceptions = {MException('correct_yasara:output_file_could_not_be_opened',...
                                    'Output file %s could not be opened',oname)};
    return;
end

chain = '';
tline = '';
while 1
    tline0 = tline;
    tline = fgetl(ifid);
    if ~ischar(tline), break, end
    if length(tline)<16 % catches too short end line
        fprintf(ofid,'%s\n',tline);
    else
        record = tline(1:6);
        if strcmp(record,'HETATM')
            record = 'ATOM  ';
        end
        if strcmp(record,'ATOM  ')
            for ks = 1:length(substitute)
                if strcmpi(tline(22),substitute(ks))
                    tline(22) = upper(newcid(ks));
                end
            end
            if ~isempty(chain)
                if ~strcmpi(chain,tline(22))
                    tline0(1:6) = 'TER   ';
                    fprintf(ofid,'%s\n',tline0(1:27));
                end
            end
            chain = tline(22);
            atname = tline(13:16);
            skip = false;
            switch atname
                case '1H5*' 
                    atname = ' H5*';
                case '2H5*' 
                    atname = 'H5**';
                case '*HO2' 
                    atname = 'HO2*';
                case ' O1P' 
                    atname = ' OP1';
                case ' O2P' 
                    atname = ' OP2';
                case '*HO3' 
                    skip = true;
                case '*HO5' 
                    skip = true;
            end
            if ~skip
                for k = 1:length(atname)
                    if atname(k) == '*'
                        atname(k) = '''';
                    end
                end
                tline(1:6) = record;
                tline(13:16) = atname;
                fprintf(ofid,'%s\n',tline);
            end
        end
    end
end

fclose(ifid);
fclose(ofid);