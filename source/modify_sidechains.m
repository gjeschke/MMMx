function [entity,modified] = modify_sidechains(entity,force,mutations)
%  entity = MODIFY_SIDECHAINS(entity,force)
%    Modify sidechains using SCWRL4
%
% INPUT:
%
% entity    input entity, if this is the only argument, sidechains with
%           missing atoms are repaired
% force     optional flag, if true, all sidechains are replaced, used
%           for complete repacking
% mutations optional list of mutations to be performed, cell of strings:
%           ctag.rtag.slc   ctag    chain tag, such as A
%                           rtag    residue tag, such as R131
%                           slc     single-letter amino acid code
%           example A.R131.H mutates residue 131 in chain A to histidine
%
% OUTPUT:
%
% entity    entity with modified sidechains
% modified  number of modified sidechains, returns with -1 if SCWRL4 is not
%           found
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

if ~exist('force','var') || isempty(force)
    force = false;
end

if ~exist('mutations','var') || isempty(mutations)
    mutations = {};
end


dospath=which('scwrl4.exe');
if isempty(dospath)
    modified = -1;
    return
end

modified = 0;
chains = fieldnames(entity);
sequence = char(zeros(1,10000));
seqpoi = 0;
to_be_replaced{1000} = ''; 
reppoi = 0;
for ch = 1:length(chains)
    chain = chains{ch};
    if isstrprop(chain(1),'upper') % chain fields start with a capital
        residues = fieldnames(entity.(chain));
        for kr = 1:length(residues) % expand over all residues
            residue = residues{kr};
            if strcmp(residue(1),'R') % these are residue fields
                slc = tlc2slc(entity.(chain).(residue).name);
                if ~isempty(slc) % this is an amino acid
                    seqpoi = seqpoi + 1;
                    complete = check_sidechain_integrity(entity.(chain).(residue));
                    if ~complete || force
                        modified = modified + 1;
                        sequence(seqpoi) = upper(slc);
                        reppoi = reppoi + 1;
                        to_be_replaced{reppoi} = [chain '.' residue];
                    else
                        sequence(seqpoi) = lower(slc);
                    end
                    for m = 1:length(mutations)
                        mutation = split(mutations{m},'.');
                        if strcmpi(mutation{1},chain) && strcmpi(mutation{2},residue) && ~isempty(mutation{3})
                            sequence(seqpoi) = upper(mutation{3});
                            if ~force && complete
                                modified = modified + 1;
                            end
                        end
                    end
                end
            end
        end
    end
end
sequence = sequence(1:seqpoi);
to_be_replaced = to_be_replaced(1:reppoi);
infile = fullfile(pwd,'SCWRL4_input.pdb');
outfile = fullfile(pwd,'SCWRL4_output.pdb');
save_options.cleaned = true;
put_pdb(entity,infile,save_options);
framefile = fullfile(pwd,'SCWRL4_frame.pdb');
save_options.cleaned = false;
save_options.hetframe = true;
put_pdb(entity,framefile,save_options);
seqfile = fullfile(pwd,'SCWRL4_sequence.seq');
fid = fopen(seqfile,'wt');
fprintf(fid,'%s\n',sequence);
fclose(fid);
comd=[dospath ' -i ' infile ' -o ' outfile ' -f ' framefile ' -s ' seqfile ' -h -t'];
[s, ~] = dos(comd);
if s~=0 
    modified = -1;
    return
end
entity2 = get_pdb(outfile);
delete(infile);
delete(outfile);
delete(framefile);
delete(seqfile);
[nat,~] = size(entity.xyz);
xyz = [entity.xyz; zeros(10000,3)];
elm = [entity.elements, zeros(1,10000),'uint8'];
occ = [entity.occupancies; zeros(10000,1,'uint8')];
for r = 1:length(to_be_replaced)
    fields = split(to_be_replaced{r},'.');
    atoms = fieldnames(entity2.(fields{1}).(fields{2}));
    for ka = 1:length(atoms) % expand over all atoms
        atom = atoms{ka};
        if isstrprop(atom(1),'upper') % these are atom fields
            % store coordinate in original entity and replace index
            index = entity2.(fields{1}).(fields{2}).(atom).tab_indices(1);
            nat = nat + 1;
            xyz(nat,:) = entity2.xyz(index,:);
            elm(nat,:) = entity2.elements(index);
            occ(nat,:) = entity2.occupancies(index);
            entity2.(fields{1}).(fields{2}).(atom).tab_indices(1) = nat;
        end
    end
    entity.(fields{1}).(fields{2}) = entity2.(fields{1}).(fields{2});
end
entity.xyz = xyz(1:nat,:);
entity.elements = elm(1:nat);
entity.occupancies = occ(1:nat);
outfile = fullfile(pwd,'SCWRL4_modified.pdb');
put_pdb(entity,outfile);
entity = get_pdb(outfile);
delete(outfile);

