function [ecoor,atomtags,code] = get_link_nt(logfid,seq,HNP_lib,anchore,anchori,ntoffset,environ)
% function [ecoor,atomtags,code] = get_link_nt(logfid,seq,HNP_lib,anchore,anchori,ntoffset,environ)
% 
% tries to find the best fitting nucleotide (RNA fragment) that links an 
% initial and a target anchor, to be used instead of mk_RNA_loop if
% a link is only one nucleotide long, both anchors must be present
%
% logfid        file handle for logfile
% seq           single letter codes from previous up to next nt, 
%               allowed nucleotides A,C,G,U
% HNP_lib       Humphris-Narayanan/Pyle pseudo-rotamer library for RNA with
%               fields 
%               .fragments          fragment library
%               .non_clash_table    table of non-clashing fragment pairs
% anchore       coordinates of the P(i+1),C4'(i+1),P(i+2) atoms of the
%               final anchor nucleotide, defaults to empty array, which is
%               interpreted as absence of a final anchor
% anchori       coordinates of the C4'(i-2), P(i-1), C4'(i-1), and O3'(i-1) 
%               atoms of the initial anchor nucleotide, defaults to empty 
%               array, which is interpreted as absence of an initial anchor
% ntoffset      nucleotide index in RNA of the preceding nucleotide,
%               defaults to 0
% environ       (ne,3) array of Cartesian coordinates of environment atoms,
%               if present and not empty, the model is tested for clashes 
%               with the environment
%
% Options:
%
% ecoor         extended coordinates, first column: nucleotide index,
%               columns 2-4: coordinates, empty if fuHNP_lib.non_clash_tableion fails
% atomtags      cell array of atom tags corresponding to ecoor
% code          fragment codes
%
% G. Jeschke, 07.05.2021

base = seq(2);

template = [anchori(3,:);anchore(1:2,:)];

transmat = eye(4);
min_rmsd = 1e6;
for kf = 1:length(HNP_lib.fragments)
    coor0 = HNP_lib.fragments(kf).(base).coor;
    atomtags = HNP_lib.fragments(kf).(base).atomtags;
    neighbours = [coor0(HNP_lib.fragments(kf).(base).assign.previous(1),:); ...
        coor0(HNP_lib.fragments(kf).(base).assign.next(2:3),:)];
    [rmsd, ~, transmat0] = superimpose_3points(template,neighbours);
    if rmsd < min_rmsd
        min_rmsd = rmsd;
        transmat = transmat0;
        coor = coor0;
        code = kf;
    end
end
[m,~] = size(coor);
coor = [coor ones(m,1)]*transmat';
coor = coor(:,1:3);

ecoor = [(1+ntoffset)*ones(m-3,1) coor(2:m-2,:)];
atomtags = atomtags(2:end-2);

fprintf(logfid,'Nucleotide inserted with neighbour r.m.s.d. of %4.1f Angstroem\n',min_rmsd);

% determine closest approach and the corresponding atom pair
if ~isempty(environ)
    pair_dist = get_all_pair_dist(coor,environ);
    min_dist = min(min(pair_dist));
    fprintf(logfid,'Minimum distance to an environment atom is %4.2f Angstroem\n',min_dist);
end

