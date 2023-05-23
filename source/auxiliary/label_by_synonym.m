function tlc = label_by_synonym(synonym)
%
% LABEL_BY_SYNOYM Returns three-letter code for a label given a synonym
%
%   tlc = LABEL_BY_SYNONYM(synonym)
%
% INPUT
% synonym       string with alternative label name, case-insensitive
%
% OUTPUT
% tlc           the three-letter code of the label in MMMx, empty if the
%               synonym is not defined
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2020: Gunnar Jeschke

% initialize empty output

tlc = '';

% Definition of known synonyms, "official" name must go first,
% command-separated, must end with a comma
synonyms{1} = 'R1A,mtsl,mtssl,';
synonyms{2} = 'R7A,br-mtsl,br-mtssl,';
synonyms{3} = 'V1A,v1,';
synonyms{4} = 'IA1,ia-proxyl,';
synonyms{5} = 'MA1,ma-proxyl,';
synonyms{6} = 'DZD,dzd,'; % placeholder only
synonyms{7} = 'DZC,dzc,'; % placeholder only
synonyms{8} = 'GDI,iag,';
synonyms{9} = 'GDM,mag,';
synonyms{10} = 'TUP,iap-4tu,iap-thiouracil,';
synonyms{11} = 'TUM,mts-4tu,mts-thiouracil,';
synonyms{12} = 'RTT,r5-tpt,';
synonyms{13} = 'R5P,r5p,'; % placeholder only
synonyms{14} = 'R3P,r3p,'; % placeholder only
synonyms{15} = 'K1H,HF-K1,';
synonyms{16} = 'NC1,cNox@Tyr,';
synonyms{17} = 'NX1,lNox@Tyr,';
synonyms{18} = 'CNR,CNC-NO,';
synonyms{19} = 'GMO,dota-gd,';
synonyms{20} = 'GTO,dtpa-gd,';
synonyms{21} = 'M8D,m8-dota-gd,';
synonyms{22} = 'GPM,gpymi-MTA,';
synonyms{23} = 'TMT,tormyshev-trityl,';
synonyms{24} = 'HCU,dHis-Cu,';
synonyms{25} = 'GBH,br-py-do3a,';
synonyms{26} = 'GBM,br-py-do3ma,';
synonyms{26} = 'HO4,HO4451,';
synonyms{27} = 'BAS,BASL,';
synonyms{28} = 'R5T,R5-TP,';

for k = 1:length(synonyms)
    if contains(upper(synonyms{k}),[',' upper(synonym) ','])...
            || strcmpi(synonyms{k}(1:3),synonym)
        tlc = synonyms{k}(1:3);
        break
    end
end
