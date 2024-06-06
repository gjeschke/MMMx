function [header,sequence] = get_fasta(inname,seqnum)
%
% function [header,sequence]=read_fasta(inname,seqnum),
%
% Reads a peptide (protein) sequence in FASTA format
%
% inname   is the filename INCLUDING extension
% seqnum   is an optional parameter that gives the number of the sequence in
%          the FASTA file, it defaults to 1
% sequence is a string that codes the sequence in single-letter format
%
% (c) G. Jeschke, 2007

    if exist('onCleanup','class')==8, % this is done only for Matlab version 
                                      % 2008a and later, where onCleanup exists
        c = onCleanup(@myCleanup);
    end;

    if nargin<2, seqnum=1; end;

    rfile=fopen(inname,'r');
    cseq=0;
    sequence=':'; % only placeholder, will be cut out at the end
    header='';

    nl=0;
    h=0;
    while h~=-1,
       line=fgetl(rfile); nl=nl+1;
       h=line;
       if h~=-1,
           if char(line(1))=='>',
               header=line;
               cseq=cseq+1; % current sequence number
               if cseq==seqnum; header=line; end;
           elseif char(line(1)) ~= ';' && cseq==seqnum, % neglect comments
               sequence=strcat(sequence,line); % append line to sequence
           end;
       end;
       if cseq>seqnum, break; end; % if we are already past the wanted sequence, stop
    end;
    if length(sequence)<2,
        sequence='';
    else
        sequence=sequence(2:length(sequence));
    end;
    fclose(rfile);

end