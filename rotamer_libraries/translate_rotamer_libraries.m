function translate_rotamer_libraries
% Translates rotamer libraries from MMM to MMMx format

R1A_tags = ':C:O:OXT:N:H2:H:CA:H4:CB:H5:H6:SG:S2:C4:H7:H8:C5:C6:H9:C7:C8:H10:H11:H12:C9:H13:H14:H15:N2:O2:C10:C11:H16:H17:H18:C12:H19:H20:H21:';
IA1_tags = ':C:O:OXT:N:H2:H:CA:H4:CB:H5:H6:SG:C4:H7:H8:C5:O2:N2:H9:C6:H10:C7:H11:H12:C8:C9:H13:H14:H15:C10:H16:H17:H18:N3:O3:C11:C12:H19:H20:H21:C13:H22:H23:H24:';
HCU_tags = ':C:O:OXT:N:H2:H:CA:H4:CB:H5:H6:CG:NE:CZ:HZ:NZ:HNZ:CE:HE:';

defs = load('labels.mat');

% load universal force field

ff = load('uff_parameters.mat');

fid = fopen('library_additions.txt','rt');
while 1
    tline = fgetl(fid);
    if ~ischar(tline) || isempty(tline)
        break 
    end
    clear rot_lib
    myline = textscan(tline,'%s');
    args=myline{1};
    rot_lib.tlc = upper(args{1});
    TLC = rot_lib.tlc;
    if strcmp(rot_lib.tlc,'K1H')
        TLC = 'RUA';
    end
    old = load(args{2});
    oname = sprintf('rotamer_library_%s.mat',args{1});
    fprintf(1,'\nCreating library %s from input %s\n',oname,args{2});
    defid = tag2id(TLC,upper(defs.label_defs.restags));
    if isempty(defid)
        fprintf(2,'Definition ID is missing for %s\n',rot_lib.tlc);
        continue
    else
        fprintf(1,'Definition ID is %i\n',defid);
    end
    rot_lib.synonyms{1} = defs.label_defs.residues(defid).short_name;
    syn1 = rot_lib.synonyms{1};
    poi = strfind(syn1,'MTSL');
    if ~isempty(poi)
        syn2 = [syn1(1:poi-1) 'MTSSL' syn1(poi+4:end)];
        rot_lib.synonyms{2} = syn2;
    end        
    if strcmpi(syn1,'M-Proxyl')
        rot_lib.synonyms{2} = 'MA-Proxyl';
    end
    iname = sprintf('%s.smi',rot_lib.tlc);
    fid_smi = fopen(iname,'rt');
    tline = fgetl(fid_smi);
    fclose(fid_smi);
    myline = textscan(tline,'%s');
    args=myline{1};
    rot_lib.SMILES = args{1};
    fprintf(1,'.SMILES: %s\n',rot_lib.SMILES);
    for r = 1:length(old.rot_lib.library)
        rot_lib.rotamers(r).coor = old.rot_lib.library(r).ecoor(:,2:4); 
        rot_lib.rotamers(r).torsion = old.rot_lib.library(r).dihedrals; 
    end
    rot_lib.elements = old.rot_lib.library(r).ecoor(:,1);
    rot_lib.populations = old.rot_lib.calibration.pop.';
    if isempty(defs.label_defs.residues(defid).spin_density)
        pos0 = old.rot_lib.usefull_atoms.midNO;
        if length(pos0) == 1
            pos = [pos0 1];
        elseif length(pos0) == 2
            pos = [pos0.' [0.5;0.5]];
        else
            fprintf(2,'Unexpected position vector size %i\n',length(pos0)); 
        end
    else
        pos = defs.label_defs.residues(defid).spin_density;
    end
    rot_lib.position = pos;
    if isfield(old.rot_lib,'attachment')
        rot_lib.attachment = old.rot_lib.attachment;
    else
        rot_lib.attachment = defs.label_defs.residues(defid).attachment;
    end
    rot_lib.side_chain = old.rot_lib.usefull_atoms.side_chain;
    if ~isfield(old.rot_lib,'forgive') || ...
            isempty(old.rot_lib.forgive)
        fprintf(1,'Assuming default forgive/attraction enhancement for label %s\n',TLC);
        rot_lib.f_factor = [0.5 1];
    else
        rot_lib.f_factor = old.rot_lib.forgive;
    end
    rot_lib.atom_tags = defs.label_defs.residues(defid).atoms; 
    mol_frame = defs.label_defs.residues(defid).frame;
    if strcmpi(TLC,'R1A')
        rot_lib.atom_tags = R1A_tags;
        mol_frame = ':N2:O2:C10:'; 
    end
    if strcmpi(TLC,'IA1')
        rot_lib.atom_tags = IA1_tags;
        mol_frame = ':N3:O3:C11:';
    end
    if strcmpi(TLC,'HCU')
        rot_lib.atom_tags = HCU_tags;
        mol_frame = ':NE:CG:CZ:';
    end
    if strcmpi(TLC,'NX1')
        mol_frame = ':NN1:ON1:C2:';
    end
    if strcmpi(TLC,'R5P')
        mol_frame = ':NO:ON:C12:';
    end
    if strcmpi(TLC,'R3P')
        mol_frame = ':NO:ON:C16:';
    end
    if strcmpi(TLC,'CNR')
        mol_frame = ':N2:O6:C9:';
    end
    if strcmpi(TLC,'GTO')
        mol_frame = ':Gd1:O4:N4:';
    end
    if strcmpi(TLC,'M8D')
        mol_frame = ':Gd1:O2:O3:';
    end
    if strcmpi(TLC,'GPM')
        mol_frame = ':Gd1:N2:N4:';
    end
    if strcmpi(TLC,'TMT')
        mol_frame = ':C19:C20:C10:';
    end
    if isempty(rot_lib.atom_tags)
        fprintf(2,'Empty atom tag list\n');
    else
        atom_tag_string = rot_lib.atom_tags;
        tag = id2tag(1,atom_tag_string);
        atom_tags = cell(1,1000);
        apoi = 0;
        while ~isempty(tag)
            apoi = apoi + 1;
            atom_tags{apoi} = tag;
            tag = id2tag(apoi+1,atom_tag_string);
        end
        rot_lib.atom_tags = atom_tags(1:apoi);
    end
    rot_lib.std_frame = old.rot_lib.stdframe;
    rot_lib.std_frame_atoms = cell(1,3);
    attach_frame = defs.label_defs.residues(defid).res_frame;
    for k = 1:3
        rot_lib.std_frame(k) = old.rot_lib.stdframe(4-k);
        if isempty(attach_frame)
            rot_lib.std_frame_atoms{k} = id2tag(rot_lib.std_frame(k),atom_tag_string);
        else
            rot_lib.std_frame_atoms{k} = id2tag(4-k,attach_frame);
        end
    end
    for k = 1:3
         tag = id2tag(k,mol_frame);
         id = tag2id(tag,atom_tag_string);
         rot_lib.mol_frame(k) = id;
    end
    rot_lib.mol_frame_atoms = mol_frame;
    if isempty(old.rot_lib.class)
        fprintf(2,'Class not defined for label %s\n',TLC);
    else
        rot_lib.class = old.rot_lib.class;
    end
    if isempty(old.rot_lib.chi_def)
        fprintf(2,'Torsion angle definition missing for label %s\n',TLC);
    else
        rot_lib.chi_def = old.rot_lib.chi_def;
    end
    if isempty(old.rot_lib.connect)
        fprintf(2,'Atom connection definition missing for label %s\n',TLC);
    else
        rot_lib.connect = old.rot_lib.connect;
    end
    if isempty(old.rot_lib.forcefield)
        fprintf(2,'Force field assignment missing for label %s\n',TLC);
    else
        rot_lib.attach_forcefield = old.rot_lib.forcefield;
    end
    if isempty(old.rot_lib.calibration.Bfactors)
        fprintf(2,'B factors missing for label %s\n',TLC);
    else
        rot_lib.B_factors = old.rot_lib.calibration.Bfactors;
    end
    if isempty(old.rot_lib.calibration.MD_program)
        fprintf(2,'Calibration method missing for label %s\n',TLC);
    else
        rot_lib.method = old.rot_lib.calibration.MD_program;
    end
    if isempty(old.rot_lib.calibration.force_field)
        fprintf(2,'Generation force field missing for label %s\n',TLC);
    else
        rot_lib.gen_forcefield = old.rot_lib.calibration.force_field;
    end
    if isempty(old.rot_lib.types)
        fprintf(2,'Atom types for force field missing for label %s\n',TLC);
    else
        rot_lib.types = old.rot_lib.types;
    end
    if isempty(old.rot_lib.solvation)
        fprintf(2,'Solvation information missing for label %s\n',TLC);
    else
        rot_lib.solvation = old.rot_lib.solvation;
    end
    if isempty(old.rot_lib.calibration.prerun)
        fprintf(2,'Prerun information missing for label %s\n',TLC);
    else
        rot_lib.prerun = old.rot_lib.calibration.prerun;
    end
    if isempty(old.rot_lib.calibration.suppress_H)
        fprintf(2,'H suppression information missing for label %s\n',TLC);
    else
        rot_lib.suppress_H = old.rot_lib.calibration.suppress_H;
    end
    if isempty(old.rot_lib.calibration.threshold)
        fprintf(2,'Acceptance threshold missing for label %s\n',TLC);
    else
        rot_lib.threshold = old.rot_lib.calibration.threshold;
    end
    if isempty(old.rot_lib.calibration.min_strain)
        fprintf(2,'Minimum strain missing for label %s\n',TLC);
    else
        rot_lib.min_strain = old.rot_lib.calibration.min_strain;
    end
    if isempty(old.rot_lib.calibration.maxdist)
        fprintf(2,'Maximum distance to backbone missing for label %s\n',TLC);
    else
        rot_lib.maxdist = old.rot_lib.calibration.maxdist;
    end
    if ~isfield(defs.label_defs.residues(defid),'color') || ...
            isempty(defs.label_defs.residues(defid).color)
        fprintf(2,'Visualization color missing for label. Assuming firebrick. %s\n',TLC);
        rot_lib.color = [178,34,34]/255;
    else
        rot_lib.color = defs.label_defs.residues(defid).color/255;
    end
    if ~isfield(defs.label_defs.residues(defid),'radius') || ...
            isempty(defs.label_defs.residues(defid).radius)
        fprintf(2,'Visualization radius missing for label. Assuming 0.64 Angstroem. %s\n',TLC);
        rot_lib.color = 0.64;
    else
        rot_lib.color = defs.label_defs.residues(defid).radius;
    end
    rot_lib.ff.LJ_r = ff.uff.LJ_r;
    rot_lib.ff.LJ_D = ff.uff.LJ_D;
    rot_lib.ff.types = ff.uff.tags;
    save(oname,'rot_lib');
end

fclose(fid);

function id=tag2id(tag,tags,codelist,delimiter)
% function id=tag2id(tag,tags,codelist,delimiter)
%
% Returns the identification code corresponding to a string tag by
% comparison with a list of possible tags given in a string with colon (:)
% separation, another separator can be selected by optional parameter
% delimiter
% if an array codelist is given the identifier is of the same type as the
% elements in codelist, the sequence in the codelist must correspond to the
% one in the tags string
% if codelist is missing, the id is an integer corresponding to the
% position of the tag in string tags
% if tag is not found in tags or if the position is larger than the number
% of elements in codelist, an empty id is given back
%
% tag       (string) tag to be found, example: 'Be'
% tags      (string) tag list, example ':H:He:Li:Be:C:N:O:F:Ne:'
% codelist  [optional] (array of ids), example [1,4,7,9,12,14,16,19,20],
%           can be empty to allow different delimiter without codelist
% delimiter [optional] (char), alternative delimiter, defaults to ':'
%
% id        id selected from array codelist or integer telling the position
%           of tag in tags, for example above:
%           if codelist is missing, id=4
%           if codelist is present, id=9
%
% G. Jeschke, 2009

id=[]; % empty default output
if nargin<4
    delimiter=':'; % colon as default delimiter
end

etag=[delimiter tag delimiter];
position=strfind(tags,etag);
if position
    id=1+sum(find(tags==delimiter)<position(1));
    if nargin>=3 && ~isempty(codelist)
        if id<=length(codelist)
            id=codelist(id);
        else
            id=[];
        end
    end
end

function tag=id2tag(id,tags,codelist,delimiter)
% function tag=id2tag(id,tags,codelist,delimiter)
%
% Returns the string tag corresponding to an identification code from
% a list of possible tags given in a string with colon (:)
% separation, another separator can be selected by optional parameter
% delimiter
% if an array codelist is given the identifier is of the same type as the
% elements in codelist, the sequence in the codelist must correspond to the
% one in the tags string
% if codelist is missing, the id is an integer corresponding to the
% position of the tag in string tags
% if the id is larger than the number of tags in the tags string an empty
% tag is given back
% if the id appears twice or more in the code list, the first tag is given
% back
%
% id        identifier, for example 7
% tags      (string) tag list, example ':H:He:Li:Be:C:N:O:F:Ne:'
% codelist  [optional] (array of ids), example [1,4,7,9,12,14,16,19,20]
% delimiter [optional] (char), alternative delimiter, defaults to ':'
%
% tag       (string) tag
%           if codelist is missing, tag for the example will be 'O'
%           if codelist is present, tag for the example will be 'Li'
%
% G. Jeschke, 2009

tag=''; % empty default output
if nargin<4
    delimiter=':'; % colon as default delimiter
end

if nargin>2
    poi=find(id==codelist);
    if ~isempty(poi)
        id=poi(end);
    else
        id=length(tags)+1;
    end
end

lookup=find(tags==delimiter);
if id<length(lookup)
    tag=tags(lookup(id)+1:lookup(id+1)-1);
end
