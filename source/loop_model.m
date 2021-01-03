function [coor,errcode,restrain,p_model,k] = loop_model(sequence, anchorN, anchorC, anchorNp, anchorCn, prot_coor, restrain, Rama_res, rescodes, n_res, options)
% [coor,errcode,restrain,p_model,k] = loop_model(sequence, anchorN, anchorC, anchorNp, anchorCn, prot_coor, restrain, Rama_res, rescodes, n_res, options)
% Generates a closed loop model conforming to residue-specific Ramachandran plots
% backbone based on: H. Sugeta, T. Miyazawa, Biopolymers, 1967, 5, 673-679.
% residue-specific Ramachandran plots from: S. Hovmöller, T. Zhou, T.
% Ohlson, Acta Cryst. D, 2002, 58, 768-776.
%
% Input:
% sequence  peptide sequence in single-letter code, must include one
%           residue before the first residue and one residue after the
%           last residue to be modelled if a closed loop, for a C-terminal
%           loop, one residue before the first until the last residue, and
%           for a non-anchored loop only the loop residues
% anchorN   backbone coordinates of N-terminal anchor residue in the order
%           N, CA, C, O (Cartesian in Å) 
% anchorC   backbone coordinates of C-terminal anchor residue in the order
%           N, CA, C, O (Cartesian in Å) 
% anchorNp  backbone coor. of the residue before the N-terminal anchor 
%           residue in the order N, CA, C, O (Cartesian in Å) 
% anchorCn  backbone coordinates of the residue after the C-terminal anchor
%           residue in the order N, CA, C, O (Cartesian in Å) 
% prot_xyz  protein coordinates for clash test
% restrain  set of restraints to be tested or enforced, loop residues,
%           format restrain(k).[type] with .[type] being
%           .secondary  0   no secondary structure
%                       1   alpha helix
%                       2   beta sheet
%                       3   inner alpha helix residue (except two residues
%                           each at the N and C terminus of the helix)
%           .aprop      propensity of alpha-helical structure (0...1)
%           .bprop      propensity of beta-strand structure (0...1)
%           .cprop      propensity of cis-peptide formation (0...1)
%           .label      label xyz coordinates relative to local frame, if
%                       empty, no distance restraints are evaluated for 
%                       this residue
%           .r_beacon   distance restraints to fixed points, vector with
%                       subfields
%                       .xyz    xyz coordinates of beacon
%                       .type   'Gaussian' or 'bounds'
%                       .par1   mean distance or lower bound, can be
%                               vectors
%                       .par2   standard deviation or upper bound,
%                               can be vectors
%           .r_intern   distance restraints to other residues in the loop,
%                       vector with subfields
%                       .site   index of other site in the loop, always
%                               smaller than k
%                       .type   'Gaussian' or 'bounds'
%                       .par1   mean distance or lower bound, can be
%                               vectors
%                       .par2   standard deviation or upper bound,
%                               can be vectors
%           .oligomer   homooligomer distance restraint, supposes that the
%                       symmetry axis is the z axis
%                       .n      number of protomers in homooligomer
%                       .type   'Gaussian' or 'bounds'
%                       .par1   mean distance or lower bound, can be
%                               vectors
%                       .par2   standard deviation or upper bound,
%                               can be vectors
%           .depth      membrane depth restraint, supposes that the
%                       membrane normal is the z axis  
%                       .site   'CA'(Calpha) or 'label
%                       .type   'Gaussian' or 'bounds'
%                       .par1   mean coordinate or lower bound, can be
%                               vectors
%                       .par2   standard deviation or upper bound,
%                               can be vectors
% Rama_res  residue_specific Ramachandran information (pairs of allowed
%           angles fi,psi)
% rescodes  MMMx-internal residue codes
% n_res     number of restraints
% options   struct with computation options
%           .min_prob   minimum restraint fulfillment probability for an
%                       accepted model, defaults to 0.5
%           .reboots    maximum number of reboots upon failure, defaults to
%                       100
%           .accstd     minimum normalized distance reductions to
%                       C-terminal anchor in second half loop, default 0.5
%           .accrad     acceptance radius for approach to N-terminal
%                       anchor, default .0
%           .clash_threshold    threshold for clashes inside loop, defaults
%                               to 2.5 Angstroem
%           .clash_threshold_lp threshold for clashes beween loop and
%                               protein, defaults to 2.5 Angstroem
%           .max_attempts       maximum attempts for generating a second
%                               half loop, defaults to 100
%           .max_res            maximum attempts for a single-residue step,
%                               defaults to 100
%
% Output:
% coor      Cartesian backbone coordinates of the loop in the order N, CA,
%           C, O; empty if not successful
% errcode   error code for unsuccessful attempts, 0 for success
%           1   second half loop could not be converged to closed loop
%           2   whole loop had a self clash
%           3   last dihedrals outside allowed Ramachandran region
%           4   whole loop clashed with protein
%           5   distance restraint violation
%           6   first half loop clashed with protein
%           7   first half loop clashed with itself
%           8   model with Ramachandran-fixed backbone violates restraints
% restrain  restraint variable with diagnostic information
%           .xyz    simulated label position
%           .r_beacon(l).p  fulfillment probability
%           .r_intern(l).p  fulfillment probability
%           .oligomer.p  fulfillment probability
%           .depth.p  fulfillment probability
% p_model   cumulative model probability from restraints
% k         residue where failure was detected
%
% G. Jeschke, 2007-2020

% set default output
errcode = 0;
p_model = 1; % default probability of an unrestrained model

if ~exist('options','var') || ~isfield(options,'min_prob')
    options.min_prob = 0.5;
end

if ~isfield(options,'reboots')
    options.reboots = 100;
end

if ~isfield(options,'accstd')
    options.accstd = 0.5;
end

if ~isfield(options,'accrad')
    options.accrad = 3.0;
end

if ~isfield(options,'clash_threshold')
    options.clash_threshold = 2.0;
end

if ~isfield(options,'clash_threshold_lp')
    options.clash_threshold_lp = 2.0;
end

if ~isfield(options,'max_attempts') % attempts to generate second half loop
    options.max_attempts = 100;
end

if ~isfield(options,'res_attempts') % attempts for single-residue step
    options.res_attempts = 100;
end

% alpha-helical region according to Hovmöller et al.
alpha_phi_LB = -89;
alpha_phi_UB = -39;
alpha_psi_LB = -66;
alpha_psi_UB = -16;
alpha_phi_psi_LB = -115; 
alpha_phi_psi_UB = -95; 

% beta_sheet region according to Hovmöller et al.
beta_phi_LB = -130;
beta_phi_UB = -105;
beta_phi_UB_proline = -80;
beta_psi_LB = 128;
beta_psi_UB = 145;

% PPII helix region
PPII_phi = -75;
PPII_psi = 160;

% cis peptide frequency according to Stewart et al., J. Mol. Biol., 1990,
% 214, 253-260

cis_Xaa_Pro = 0.065;
cis_Xaa_non_Pro = 0.0003;

dextended = 3.8; % extended mean Calpha-Calpha distance
dhelix = 1.5; % helical (contracted) mean Calpha-Calpha distance
dmean = (dextended + dhelix)/2; % mean Calpha-Calpha distance step

if ~isempty(anchorC)
    sequence = sequence(1:end-1);
end

% Default bond lengths
b1=1.416; % N-C_alpha
b1vec=[b1,0,0];
b2=1.529; % C_alpha-C
b2vec=[b2,0,0];
b3=1.320; % C-N
b3vec=[b3,0,0];

% Default bond angles
fi1=110*pi/180; % N-C_alpha-C
sfi1 = sin(fi1); cfi1 = cos(fi1);
fi2=114.2*pi/180; % C_alpha-C-N
sfi2 = sin(fi2); cfi2 = cos(fi2);
fi3=121*pi/180; % C-N-C_alpha
sfi3 = sin(fi3); cfi3 = cos(fi3);
omega_trans = pi;
sot = sin(omega_trans); cot = cos(omega_trans);
omega_cis = 0;
soc = sin(omega_cis); coc = cos(omega_cis);

% process propensity restraints, if any
if isfield(restrain,'aprop')
    for k = 1:length(restrain)
        if restrain(k).secondary == 0 && ~isempty(restrain(k).aprop)
            if rand <= restrain(k).aprop
                restrain(k).secondary = 1;
            end
        end
    end
    % determine 'inside helix' residues
    count = 0;
    for k = 1:length(restrain)
        if restrain(k).secondary == 1
            count = count + 1;
        end
        if restrain(k).secondary ~= 1 && restrain(k).secondary ~= 3
            count = 0;
        end
        if count > 2
            restrain(k).secondary = 3;
        end
    end
    count = 0;
    for k = length(restrain):-1:1
        if restrain(k).secondary == 3
            count = count + 1;
        end
        if restrain(k).secondary ~= 1 && restrain(k).secondary ~= 3
            count = 0;
        end
        if count < 3 && count > 0
            restrain(k).secondary = 1;
        end
    end
end

if isfield(restrain,'bprop')
    for k = 1:length(restrain)
        if restrain(k).secondary == 0 && ~isempty(restrain(k).bprop)
            if rand <= restrain(k).bprop
                restrain(k).secondary = 2;
            end
        end
    end
end

if isfield(restrain,'pprop')
    for k = 1:length(restrain)
        if restrain(k).secondary == 0 && ~isempty(restrain(k).pprop)
            if rand <= restrain(k).pprop
                restrain(k).secondary = 4;
            end
        end
    end
end

if ~isempty(anchorC)
    ngap = length(sequence)-2;
    backbone = zeros(4*(ngap+3),3);
    backbone(end-7:end-4,:) = anchorC;
    backbone(end-3:end,:) = anchorCn;
else
    backbone = zeros(4*(length(sequence)+1),3);
    ngap = length(sequence)-1;
end

phivec = zeros(1,ngap);
psivec = zeros(1,ngap);

nout = floor(ngap/2);

rejections = zeros(1,ngap+2);
drvec = rejections;
dRvec = rejections;

% initial coordinates, if N-terminal anchor is not defined
acoor=[0;0;0];
A=eye(3);
k = 2;

free_standing = true;
if ~isempty(anchorN)
    free_standing = false;
    backbone(1:4,:) = anchorN;

    % initial coordinates with N-terminal anchor, psi of N terminal anchor not yet defined
    C_back = anchorNp(3,:);
    N = anchorN(1,:);
    CA = anchorN(2,:);
    C = anchorN(3,:);
    phi = dihedral_fast(C_back,N,CA,C);

    poi = round(rand*Rama_res.me{rescodes(k)}+0.5);
    ephi=Rama_res.ephi{rescodes(k)}(poi);
    epsi=Rama_res.epsi{rescodes(k)}(poi);

    [~,poi] = min(abs(ephi-phi)); % find close match of N-terminal anchor phi
    psi = epsi(poi); % set a corresponding psi value

    % local frame of Calpha of the N-terminal anchor residue
    x = anchorN(3,:)-anchorN(2,:);
    bondx = norm(x);
    x = x/norm(x);
    yp = anchorN(2,:)-anchorN(1,:);
    yp = yp/norm(yp);
    z = cross_rowvec(x,yp);
    z = z/norm(z);
    y = cross_rowvec(z,x);
    y = y/norm(y);
    A = [x;y;z];
    % transformation matrix into that local frame
    A = A';
    % coordinates of Calpha of the N-terminal anchor residue
    acoor = anchorN(2,:);
    acoor = acoor';
    acoor = acoor + A*[bondx;0;0]; % coordinates of C of the N-terminal anchor residue

    ctau=cos(psi);
    stau=sin(psi);
    A23=[-cfi2,-sfi2,0;sfi2*ctau,-cfi2*ctau,-stau;sfi2*stau,-cfi2*stau,ctau];
    A=A*A23;
    acoor=acoor+A*b3vec'; % coordinates of N of first gap residue

    % decide whether cis or trans peptide
    so = sot;
    co = cot;
    dice = rand;
    if upper(sequence(2))~='P'
        if dice < cis_Xaa_non_Pro
            so = soc; co = coc;
        end
    else
        if dice < cis_Xaa_Pro
            so = soc; co = coc;
        end
    end
    A31=[-cfi3,-sfi3,0;sfi3*co,-cfi3*co,-so;sfi3*so,-cfi3*so,co];
    A=A*A31; % local frame at N of first gap residue
    k = 2;
end

kend1 = nout+1;
if isempty(anchorC)
    kend1 = length(sequence);
    sequence = strcat(sequence,'.');
    if free_standing
        kend1 = kend1 + 1;
        sequence = strcat('.',sequence);
        rescodes = [0 rescodes];
    end
end

reboot.k = k - 1;
reboot.counter = reboots;
reboot.A = A;
reboot.acoor = acoor;
reboot.p_model = p_model;
reboot.tested = 0;

reboot_thresh = options.min_prob^(1/n_res);
tested_restraints = 0;

% first half loop generation, index k starts with k = 1 at existing anchor 
% for C-terminal loop, this generates the whole loop
while k <= kend1
    updated = false;
    if restrain(k-1).secondary == 3
        [rphi,rpsi] = get_phi_psi_in_helix(sequence(k));
    elseif restrain(k-1).secondary == 4
        phi = PPII_phi + 5*randn;
        rphi = pi*phi/180;
        psi = PPII_psi + 5*randn;
        rpsi = pi*psi/180;
    else
        poi = round(rand*Rama_res.me{rescodes(k)}+0.5);
        rphi=Rama_res.ephi{rescodes(k)}(poi);
        rpsi=Rama_res.epsi{rescodes(k)}(poi);
    end
    phivec(k) = 180*rphi/pi;
    psivec(k) = 180*rpsi/pi;
    so = sot;
    co = cot;
    dice = rand;
    if upper(sequence(k+1))~='P'
        if dice < cis_Xaa_non_Pro
            so = soc; co = coc;
        end
    else
        if dice < cis_Xaa_Pro
            so = soc; co = coc;
        end
    end
    if restrain(k-1).cis
        so = soc;
        co = coc;
    end
    if restrain(k-1).secondary
        att_sec = 1;
        if restrain(k-1).secondary == 1 % alpha-helix
            while phivec(k) < alpha_phi_LB || phivec(k) > alpha_phi_UB || ...
                  psivec(k) < alpha_psi_LB || psivec(k) > alpha_psi_UB || ...
                  phivec(k) + psivec(k) < alpha_phi_psi_LB || ...
                  phivec(k) + psivec(k) > alpha_phi_psi_UB
                    poi = round(rand*Rama_res.me{rescodes(k)}+0.5);
                    rphi=Rama_res.ephi{rescodes(k)}(poi);
                    rpsi=Rama_res.epsi{rescodes(k)}(poi);
                    phivec(k) = 180*rphi/pi;
                    psivec(k) = 180*rpsi/pi;
            end
        elseif restrain(k-1).secondary == 2 % beta-sheet
            if upper(sequence(k+1))~='P'
                while phivec(k) < beta_phi_LB || phivec(k) > beta_phi_UB || ...
                      psivec(k) < beta_psi_LB || psivec(k) > beta_psi_UB
                        poi = round(rand*Rama_res.me{rescodes(k)}+0.5);
                        rphi=Rama_res.ephi{rescodes(k)}(poi);
                        rpsi=Rama_res.epsi{rescodes(k)}(poi);
                        phivec(k) = 180*rphi/pi;
                        psivec(k) = 180*rpsi/pi;
                        att_sec = att_sec + 1;
                end
            else % special handling for proline case
                while phivec(k) < beta_phi_LB || phivec(k) > beta_phi_UB_proline || ...
                      psivec(k) < beta_psi_LB || psivec(k) > beta_psi_UB
                        poi = round(rand*Rama_res.me{rescodes(k)}+0.5);
                        rphi=Rama_res.ephi{rescodes(k)}(poi);
                        rpsi=Rama_res.epsi{rescodes(k)}(poi);
                        phivec(k) = 180*rphi/pi;
                        psivec(k) = 180*rpsi/pi;
                        att_sec = att_sec + 1;
                end
            end
        end
    end
    backbone(4*k-3,:)=acoor'; % coordinates of N
    acoor=acoor+A*b1vec'; % coordinates of C_alpha
    backbone(4*k-2,:)=acoor';
    ctau=cos(rphi);
    stau=sin(rphi);
    A12=[-cfi1,-sfi1,0;sfi1*ctau,-cfi1*ctau,-stau;sfi1*stau,-cfi1*stau,ctau];
    A=A*A12;
    acoor=acoor+A*b2vec'; % coordinates of C
    backbone(4*k-1,:)=acoor';
    if ~isempty(restrain(k-1).label) 
        % make spin label coordinate
        x= backbone(4*k-3,:) - backbone(4*k-2,:); % x axis is along C_alpha-N bond
        x=x/norm(x);    % unit vector along x
        yp=backbone(4*k-1,:) - backbone(4*k-2,:); % y axis is in the plane spanned by x axis and C-Ca bond
        yp=yp/norm(yp);
        z=cross_rowvec(x,yp); % z axis is perpendicular on xy plane
        z=z/norm(z);
        y=cross_rowvec(z,x); % real (corrected) y axis 
        dircos=[x;y;z];
        Rp=dircos; % rotation matrix for conversion to standard frame
        restrain(k-1).xyz = restrain(k-1).label*Rp + backbone(4*k-2,:);
        p_beacon = 1;
        for kr = 1:length(restrain(k-1).r_beacon)
            r = norm(restrain(k-1).xyz-restrain(k-1).r_beacon(kr).xyz); 
            switch restrain(k-1).r_beacon(kr).type
                case 'Gaussian'
                    restrain(k-1).r_beacon(kr).p = ...
                        exp(-((r-restrain(k-1).r_beacon(kr).par1)/(sqrt(2)*restrain(k-1).r_beacon(kr).par2))^2);
                    p_beacon = p_beacon * restrain(k-1).r_beacon(kr).p;
                    tested_restraints = tested_restraints + 1;
                    updated = true;
                case 'bounds'
                    if r < restrain(k-1).r_beacon(kr).par1 || ...
                            r > restrain(k-1).r_beacon(kr).par2 % bounds violated
                       restrain(k-1).r_beacon(kr).p = 0;
                       p_beacon = 0;
                    end                        
                otherwise
                    error('MMMx:loop_model:unknownRestraintType','Restraint type %s not known',restrain(k-1).r_beacon(kr).type);
            end
        end
        p_model = p_beacon*p_model;
        if p_model < min_prob
            if reboot.counter > 0
                reboot.counter = reboot.counter - 1;
                A = reboot.A;
                acoor = reboot.acoor;
                k = reboot.k + 1;
                tested_restraints = reboot.tested;
                p_model = reboot.p_model;
                continue
            end
            errcode = 5;
            coor = [];
            return
        end
        p_intern = 1;
        for kr = 1:length(restrain(k-1).r_intern)
            site = restrain(k-1).r_intern(kr).site;
            r = norm(restrain(k-1).xyz-restrain(site).xyz); 
            switch restrain(k-1).r_intern(kr).type
                case 'Gaussian'
                    restrain(k-1).r_intern(kr).p = ...
                        exp(-((r-restrain(k-1).r_intern(kr).par1)/(sqrt(2)*restrain(k-1).r_intern(kr).par2))^2);
                    p_intern = p_intern * restrain(k-1).r_intern(kr).p;                 
                    tested_restraints = tested_restraints + 1;
                    updated = true;
                case 'bounds'
                    if r < restrain(k-1).r_intern(kr).par1 || ...
                            r > restrain(k-1).r_intern(kr).par2 % bounds violated
                        restrain(k-1).r_intern(kr).p = 0;
                        p_intern = 0;
                    end                        
                otherwise
                    error('MMMx:loop_model:unknownRestraintType','Restraint type %s not known',restrain(k-1).r_beacon(kr).type);
            end
        
        end
        p_model = p_intern*p_model;
        if p_model < options.min_prob
            if reboot.counter > 0
                reboot.counter = reboot.counter - 1;
                A = reboot.A;
                acoor = reboot.acoor;
                k = reboot.k + 1;
                tested_restraints = reboot.tested;
                p_model = reboot.p_model;
                continue
            end
            errcode = 5;
            coor = [];
            return
        end
        p_oligomer = 1;
        for kr = 1:length(restrain(k-1).oligomer)
            r = 2*sqrt(sum(restrain(k-1).xyz(1:2).^2))*sin(pi/restrain(k-1).oligomer(kr).n); 
            switch restrain(k-1).oligomer.type
                case 'Gaussian'
                    restrain(k-1).oligomer(kr).p = ...
                        exp(-((r-restrain(k-1).oligomer(kr).par1)/(sqrt(2)*restrain(k-1).oligomer(kr).par2))^2);
                    p_oligomer = p_oligomer * restrain(k-1).oligomer(kr).p;
                    tested_restraints = tested_restraints + 1;
                    updated = true;                     
                case 'bounds'
                    if r < restrain(k-1).oligomer(kr).par1 || ...
                            r > restrain(k-1).oligomer.par2 % bounds violated
                        restrain(k-1).oligomer(kr).p = 0;
                        p_oligomer = 0;
                    end                        
                otherwise
                    error('MMMx:loop_model_reverse:unknownRestraintType','Restraint type %s not known',restrain(k).oligomer(kr).type);
            end
        end
        p_model = p_oligomer*p_model;
        if p_model < min_prob
            if reboot.counter > 0
                reboot.counter = reboot.counter - 1;
                A = reboot.A;
                acoor = reboot.acoor;
                k = reboot.k + 1;
                tested_restraints = reboot.tested;
                p_model = reboot.p_model;
                continue
            end
            errcode = 5;
            coor = [];
            return
        end
        p_depth = 1;
        for kr = 1:length(restrain(k-1).depth)
            switch restrain(k-1).depth(kr).site
                case 'CA'
                    z = abs(backbone(4*k-2,3));
                case 'label'
                    z = abs(restrain(k-1).xyz(3));
                otherwise
                    error('MMMx:loop_model:unknownRestraintModifier','Depth restraint site %s not known',restrain(k).depth(kr).site);
            end
            switch restrain(k-1).depth.type
                case 'Gaussian'
                    restrain(k-1).depth(kr).p = ...
                        exp(-((z-restrain(k-1).depth(kr).par1)/(sqrt(2)*restrain(k-1).depth(kr).par2))^2);
                    p_depth = p_depth * restrain(k-1).depth(kr).p;                 
                    tested_restraints = tested_restraints + 1;
                    updated = true;
                case 'bounds'
                    if z >= restrain(k-1).depth(kr).par1 && ...
                            z <= restrain(k-1).depth.par2 % bounds not violated
                        restrain(k-1).depth(kr).p = 1;
                        p_depth = 1;
                    else
                        restrain(k-1).depth(kr).p = 1;
                        p_depth = 0;
                    end                        
                otherwise
                    error('MMMx:loop_model:unknownRestraintType','Restraint type %s not known',restrain(k).depth(kr).type);
            end
        end
        p_model = p_depth*p_model;
        if p_model < min_prob
            if reboot.counter > 0
                reboot.counter = reboot.counter - 1;
                A = reboot.A;
                acoor = reboot.acoor;
                k = reboot.k + 1;
                tested_restraints = reboot.tested;
                p_model = reboot.p_model;
                continue
            end
            errcode = 5;
            coor = [];
            return
        end
    end
    ctau=cos(rpsi);
    stau=sin(rpsi);
    A23=[-cfi2,-sfi2,0;sfi2*ctau,-cfi2*ctau,-stau;sfi2*stau,-cfi2*stau,ctau];
    A=A*A23;
    acoor=acoor+A*b3vec'; % coordinates of next N
    backbone(4*k+1,:)=acoor';
    A31=[-cfi3,-sfi3,0;sfi3*co,-cfi3*co,-so;sfi3*so,-cfi3*so,co];
    A=A*A31;    
    if p_model^(1/tested_restraints) > reboot_thresh && updated % set reboot point, if model is 'better than expected' at this poiunt
        reboot.k = k;
        reboot.counter = reboots;
        reboot.A = A;
        reboot.acoor = acoor;
        reboot.p_model = p_model;
        reboot.tested = tested_restraints;
    end
    k = k + 1;
end

% add backbone O atoms
for k = 2:kend1
    backbone = add_O(k,backbone);
end

if ~isempty(prot_coor)
    % check for clashes with protein
    pair_dist = get_all_pair_dist(backbone(9:4*(kend1-1)+3,:),prot_coor);
    min_dist = min(min(pair_dist));
    if min_dist < clash_threshold_lp
        errcode = 6;
        coor = [];
        return
    end
end

% self-clash test
k = 1;
clash = false;
while ~clash && k <= kend1-2
    min_dist = get_min_pair_dist(k,backbone(1:4*(kend1-1)+4,:));
    if min_dist < clash_threshold
        clash = true;
    end
    k = k + 1;
end

if clash
    errcode = 7;
    coor = [];
    return
end

if isempty(anchorC) % C-terminal loop is complete at this point
    coor = backbone(5:end,:); % exclude N-terminal anchor
    return;
end

k = nout+1;
A00 = A;
acoor00 = acoor;
failed = 0;
p_model_half = p_model;

if p_model^(1/tested_restraints) > reboot_thresh
    reboot.k = k;
    reboot.counter = reboots;
    reboot.A = A;
    reboot.acoor = acoor;
    reboot.p_model = p_model;
    reboot.tested = tested_restraints;
else
    reboot.counter = 0;
end
        
% second half loop generation
while k<ngap+1 && failed < options.max_attempts
    updated = false;
    Rvec = anchorC(1,:) - acoor';
    R = norm(Rvec);
    dR = R/(ngap - k + 1); % required mean distance stept towards C anchor
    k = k + 1;
    unaccepted = true;
    attempts = 0;
    acoor0 = acoor;
    A0 = A;
    detvec = zeros(1,options.res_attempts);
    while unaccepted && attempts < options.res_attempts
        attempts = attempts + 1;
        if restrain(k-1).secondary == 3
            [rphi,rpsi] = get_phi_psi_in_helix(sequence(k));
        elseif restrain(k-1).secondary == 4
            phi = PPII_phi + 5*randn;
            rphi = pi*phi/180;
            psi = PPII_psi + 5*randn;
            rpsi = pi*psi/180;
        else
            poi = round(rand*Rama_res.me{rescodes(k)}+0.5);
            rphi=Rama_res.ephi{rescodes(k)}(poi);
            rpsi=Rama_res.epsi{rescodes(k)}(poi);
        end
        phivec(k) = 180*rphi/pi;
        psivec(k) = 180*rpsi/pi;
        so = sot;
        co = cot;
        dice = rand;
        if upper(sequence(k+1))~='P'
            if dice < cis_Xaa_non_Pro
                so = soc; co = coc;
            end
        else
            if dice < cis_Xaa_Pro
                so = soc; co = coc;
            end
        end
        if restrain(k-1).cis
            so = soc;
            co = coc;
        end
        if restrain(k-1).secondary
            att_sec = 1;
            if restrain(k-1).secondary == 1 % alpha-helix
                while phivec(k) < alpha_phi_LB || phivec(k) > alpha_phi_UB || ...
                      psivec(k) < alpha_psi_LB || psivec(k) > alpha_psi_UB || ...
                      phivec(k) + psivec(k) < alpha_phi_psi_LB || ...
                      phivec(k) + psivec(k) > alpha_phi_psi_UB
                        poi = round(rand*Rama_res.me{rescodes(k)}+0.5);
                        rphi=Rama_res.ephi{rescodes(k)}(poi);
                        rpsi=Rama_res.epsi{rescodes(k)}(poi);
                        phivec(k) = 180*rphi/pi;
                        psivec(k) = 180*rpsi/pi;
                end
            elseif restrain(k-1).secondary == 2 % beta-sheet
                if upper(sequence(k+1))~='P'
                    while phivec(k) < beta_phi_LB || phivec(k) > beta_phi_UB || ...
                          psivec(k) < beta_psi_LB || psivec(k) > beta_psi_UB
                            poi = round(rand*Rama_res.me{rescodes(k)}+0.5);
                            rphi=Rama_res.ephi{rescodes(k)}(poi);
                            rpsi=Rama_res.epsi{rescodes(k)}(poi);
                            phivec(k) = 180*rphi/pi;
                            psivec(k) = 180*rpsi/pi;
                            att_sec = att_sec + 1;
                    end
                else % special handling for proline case
                    while phivec(k) < beta_phi_LB || phivec(k) > beta_phi_UB_proline || ...
                          psivec(k) < beta_psi_LB || psivec(k) > beta_psi_UB
                            poi = round(rand*Rama_res.me{rescodes(k)}+0.5);
                            rphi=Rama_res.ephi{rescodes(k)}(poi);
                            rpsi=Rama_res.epsi{rescodes(k)}(poi);
                            phivec(k) = 180*rphi/pi;
                            psivec(k) = 180*rpsi/pi;
                            att_sec = att_sec + 1;
                    end
                end
            end
        end
        backbone(4*k-3,:)=acoor'; % coordinates of N
        acoor=acoor+A*b1vec'; % coordinates of C_alpha
        backbone(4*k-2,:)=acoor';
        ctau=cos(rphi);
        stau=sin(rphi);
        A12=[-cfi1,-sfi1,0;sfi1*ctau,-cfi1*ctau,-stau;sfi1*stau,-cfi1*stau,ctau];
        A=A*A12;
        acoor=acoor+A*b2vec'; % coordinates of C
        backbone(4*k-1,:)=acoor';
        ctau=cos(rpsi);
        stau=sin(rpsi);
        A23=[-cfi2,-sfi2,0;sfi2*ctau,-cfi2*ctau,-stau;sfi2*stau,-cfi2*stau,ctau];
        A=A*A23;
        acoor=acoor+A*b3vec'; % coordinates of next N
        A31=[-cfi3,-sfi3,0;sfi3*co,-cfi3*co,-so;sfi3*so,-cfi3*so,co];
        A=A*A31;    
        rvec = acoor' - acoor0';
        dr = sum(rvec.*Rvec)/R; % distance change towards C anchor
        detm = (dr - dR)/dmean;
        detvec(attempts) = detm;
        if abs(detm) < options.accstd
            unaccepted = false;
        elseif attempts < options.res_attempts
            acoor = acoor0;
            A = A0;
        end
    end
    rejections(k) = attempts-1;
    drvec(k) = dr;
    dRvec(k) = dR;
    if k >= ngap+1 && norm(acoor'-anchorC(1,:)) > options.accrad
        A = A00;
        acoor = acoor00;
        k = nout + 1;
        failed = failed +1;
        p_model = p_model_half;
    end
    if ~isempty(restrain(k-1).label) 
        % make spin label coordinate
        x= backbone(4*k-3,:) - backbone(4*k-2,:); % x axis is along C_alpha-N bond
        x=x/norm(x);    % unit vector along x
        yp=backbone(4*k-1,:) - backbone(4*k-2,:); % y axis is in the plane spanned by x axis and C-Ca bond
        yp=yp/norm(yp);
        z=cross_rowvec(x,yp); % z axis is perpendicular on xy plane
        z=z/norm(z);
        y=cross_rowvec(z,x); % real (corrected) y axis 
        dircos=[x;y;z];
        Rp=dircos; % rotation matrix for conversion to standard frame
        restrain(k-1).xyz = restrain(k-1).label*Rp + backbone(4*k-2,:);
        p_beacon = 1;
        for kr = 1:length(restrain(k-1).r_beacon)
            r = norm(restrain(k-1).xyz-restrain(k-1).r_beacon(kr).xyz); 
            switch restrain(k-1).r_beacon(kr).type
                case 'Gaussian'
                    restrain(k-1).r_beacon(kr).p = ...
                        exp(-((r-restrain(k-1).r_beacon(kr).par1)/(sqrt(2)*restrain(k-1).r_beacon(kr).par2))^2);
                    p_beacon = p_beacon * restrain(k-1).r_beacon(kr).p;                 
                    tested_restraints = tested_restraints + 1;
                    updated = true;
                case 'bounds'
                    if r < restrain(k-1).r_beacon(kr).par1 || ...
                            r > restrain(k-1).r_beacon(kr).par2 % bounds violated
                        restrain(k-1).r_beacon(kr).p = 0;
                        p_beacon = 0;
                    end                       
                otherwise
                    error('MMMx:loop_model:unknownRestraintType','Restraint type %s not known',restrain(k-1).r_beacon(kr).type);
            end
        end
        p_model = p_beacon*p_model;
        if p_model < min_prob
            if reboot.counter > 0
                reboot.counter = reboot.counter - 1;
                A = reboot.A;
                acoor = reboot.acoor;
                k = reboot.k;
                tested_restraints = reboot.tested;
                p_model = reboot.p_model;
                continue
            end
            errcode = 5;
            coor = [];
%                 fprintf(2,'Model rejected at residue %i by too low restraint fulfillment (beacon).\n',k-1);
            return
        end
        p_intern = 1;
        for kr = 1:length(restrain(k-1).r_intern)
            site = restrain(k-1).r_intern(kr).site;
            r = norm(restrain(k-1).xyz-restrain(site).xyz); 
            switch restrain(k-1).r_intern(kr).type
                case 'Gaussian'
                    p_intern = p_intern * ...
                     exp(-((r-restrain(k-1).r_intern(kr).par1)/(sqrt(2)*restrain(k-1).r_intern(kr).par2))^2);
                    tested_restraints = tested_restraints + 1;
                    updated = true;
               case 'bounds'
                    if r < restrain(k-1).r_intern(kr).par1 || ...
                            r > restrain(k-1).r_intern(kr).par2 % bounds violated
                        restrain(k-1).r_intern(kr).p = 0;
                        p_intern = 0;
                    end                       
                otherwise
                    error('MMMx:loop_model:unknownRestraintType','Restraint type %s not known',restrain(k-1).r_beacon(kr).type);
            end
        end
        p_model = p_intern*p_model;
        if p_model < options.min_prob
            if reboot.counter > 0
                reboot.counter = reboot.counter - 1;
                A = reboot.A;
                acoor = reboot.acoor;
                k = reboot.k;
                tested_restraints = reboot.tested;
                p_model = reboot.p_model;
                continue
            end
            errcode = 5;
            coor = [];
            return
        end
        p_oligomer = 1;
        for kr = 1:length(restrain(k-1).oligomer)
            r = 2*sqrt(sum(restrain(k-1).xyz(1:2).^2))*sin(pi/restrain(k-1).oligomer(kr).n); 
            switch restrain(k-1).oligomer.type
                case 'Gaussian'
                    restrain(k-1).oligomer(kr).p = ...
                        exp(-((r-restrain(k-1).oligomer(kr).par1)/(sqrt(2)*restrain(k-1).oligomer(kr).par2))^2);
                    p_oligomer = p_oligomer * restrain(k-1).oligomer(kr).p;
                    tested_restraints = tested_restraints + 1;
                    updated = true;                    
                case 'bounds'
                    if r < restrain(k-1).oligomer(kr).par1 || ...
                            r > restrain(k-1).oligomer.par2 % bounds violated
                        restrain(k-1).oligomer(kr).p = 0;
                        p_oligomer = 0;
                    end                        
                otherwise
                    error('MMMx:loop_model:unknownRestraintType','Restraint type %s not known',restrain(k-1).oligomer(kr).type);
            end
        end
        p_model = p_oligomer*p_model;
        if p_model < options.min_prob
            if reboot.counter > 0
                reboot.counter = reboot.counter - 1;
                A = reboot.A;
                acoor = reboot.acoor;
                k = reboot.k;
                tested_restraints = reboot.tested;
                p_model = reboot.p_model;
                continue
            end
            errcode = 5;
            coor = [];
            return
        end
        p_depth = 1;
        for kr = 1:length(restrain(k-1).depth)
            switch restrain(k-1).depth(kr).site
                case 'CA'
                    z = abs(backbone(4*k-2,3));
                case 'label'
                    z = abs(restrain(k-1).xyz(3));
                otherwise
                    error('MMMx:loop_model:unknownRestraintModifier','Depth restraint site %s not known',restrain(k-1).depth(kr).site);
            end
            switch restrain(k-1).depth.type
                case 'Gaussian'
                    restrain(k-1).depth(kr).p = ...
                        exp(-((z-restrain(k-1).depth(kr).par1)/(sqrt(2)*restrain(k-1).depth(kr).par2))^2);
                    p_depth = p_depth * restrain(k-1).depth(kr).p;                 
                    tested_restraints = tested_restraints + 1;
                    updated = true;
                case 'bounds'
                    if z >= restrain(k-1).depth(kr).par1 && ...
                            z <= restrain(k-1).depth.par2 % bounds not violated
                        restrain(k-1).depth(kr).p = 1;
                        p_depth = 1;
                    else
                        p_depth = 0;
                    end                        
                otherwise
                    error('MMMx:loop_model:unknownRestraintType','Restraint type %s not known',restrain(k-1).depth(kr).type);
            end
        end
        p_model = p_depth*p_model;
        if p_model < options.min_prob
            if reboot.counter > 0
                reboot.counter = reboot.counter - 1;
                A = reboot.A;
                acoor = reboot.acoor;
                k = reboot.k;
                tested_restraints = reboot.tested;
                p_model = reboot.p_model;
                continue
            end
            errcode = 5;
            coor = [];
            return
        end
    end
    if p_model^(1/tested_restraints) > reboot_thresh && updated % set reboot point, if model is 'better than expected' at this point
        reboot.k = k;
        reboot.counter = reboots;
        reboot.A = A;
        reboot.acoor = acoor;
        reboot.p_model = p_model;
        reboot.tested = tested_restraints;
    end
end

kres = k;
if failed >= options.max_attempts
    errcode = 1;
    coor = [];
    return
else
end
% redetermine dihedrals
rphivec = zeros(1,ngap+2);
rpsivec = rphivec;
romvec = rphivec;
for k = 2:ngap+1
    [phi,psi,omega] = dihedrals(k,backbone);
    rphivec(k) = 180*phi/pi;
    rpsivec(k) = 180*psi/pi;
    romvec(k) = 180*omega/pi;
end
backbone(4*(ngap+2)-3,:) = acoor';
acoor=acoor+A*b1vec'; % coordinates of C_alpha
backbone(4*(ngap+2)-2,:) = acoor';

% stretch backbone to make C-terminal anchor CA overlap
R0 = anchorC(2,:) - acoor'; % difference vector
steps = 3*ngap + 1; 
dR0 = R0/steps;
for k = 1:steps
    res = 2 + floor(k/3); % number of residue to be corrected, starts with first gap residue 2
    poi = 4*(res-1) + 1 + mod(k,3);
    backbone(poi,:) = backbone(poi,:) + k*dR0;
end

for k = 2:ngap+1
    backbone = add_O(k,backbone);
end

% self-clash test
k = 1;
clash = false;
while ~clash && k <= ngap
    min_dist = get_min_pair_dist(k,backbone);
    if min_dist < options.clash_threshold
        clash = true;
    end
    k = k + 1;
end

k = kres;

if clash
    errcode = 2;
    coor = [];
    return
end
coor = backbone(5:end,:); % exclude N-terminal anchor

allowed = test_Rama(ngap+2,backbone,Rama_res,sequence);
if ~allowed
    [backbone_f,success] = try_fix_rama(ngap+2,backbone,Rama_res,sequence,+1);
    if ~success
        [backbone_f,success] = try_fix_rama(ngap+2,backbone,Rama_res,sequence,-1);
    end
    if ~success
        errcode = 3;
        return
    else
        backbone = backbone_f;
        errcode = test_constraints(backbone,restrain,options,ngap+1);
        if errcode
            errcode = 8;
            coor = [];
            return
        else
            errcode = -1;
        end
    end
end

pair_dist = get_all_pair_dist(backbone(9:4*ngap,:),prot_coor);
min_dist = min(min(pair_dist));
if min_dist < options.clash_threshold_lp
    if errcode == -1
        errcode = -4;
    else
        errcode = 4;
    end
    return
else
end

function [phi,psi,omega] = dihedrals(k,backbone)

[m,~] = size(backbone);
CA_back = backbone(4*k-6,:);
C_back = backbone(4*k-5,:);
N = backbone(4*k-3,:);
CA = backbone(4*k-2,:);
C = backbone(4*k-1,:);
omega = dihedral_fast(CA_back,C_back,N,CA);
phi = dihedral_fast(C_back,N,CA,C);
if 4*k+1 > m
    psi = [];
else
    N_next = backbone(4*k+1,:);
    psi = dihedral_fast(N,CA,C,N_next);
end

function backbone = add_O(k,backbone)
% carbonyl O position in standard frame (CA-N is x axis, CA-C in xy plane,
% CA is origin)

CA = backbone(4*k-2,:);
C = backbone(4*k-1,:);
Nn = backbone(4*k+1,:);

rO = 1.22;

v1 = CA-C;
v2 = Nn-C;
v1n = v1/norm(v1);
v2n = v2/norm(v2);
z = cross_rowvec(v1n,v2n);
z = z/norm(z);
Rp = rotmatn(2*pi/3,z);
backbone(4*k,:) = C + rO*v1n*Rp;

function transmat = rotmatn(theta,n)

transmat = zeros(3);
c=cos(theta);
s=sin(theta(1));
t=1-c;
n=n/norm(n);
nx=n(1);
ny=n(2);
nz=n(3);
transmat(1,1)=t*nx^2+c;
transmat(1,2)=t*nx*ny-s*nz;
transmat(1,3)=t*nx*nz+s*ny;
transmat(2,1)=t*nx*ny+s*nz;
transmat(2,2)=t*ny^2+c;
transmat(2,3)=t*ny*nz-s*nx;
transmat(3,1)=t*nx*nz-s*ny;
transmat(3,2)=t*ny*nz+s*nx;
transmat(3,3)=t*nz^2+c;

function min_dist = get_min_pair_dist(k,backbone)

a = backbone(4*k-3:4*k,:);
b = backbone(4*k+5:end,:);

pair_dist = get_all_pair_dist(a,b);
min_dist = min(min(pair_dist));

function pair_dist = get_all_pair_dist(a,b)

[m1,~] = size(a); % get sizes of the coordinates arrays
[m2,~] = size(b);

a2 = repmat(sum(a.^2,2),1,m2);
b2 = repmat(sum(b.^2,2),1,m1).';
pair_dist = sqrt(abs(a2 + b2 - 2*a*b.'));

function [phi,psi] = get_phi_psi_in_helix(res)
% secondary structure restraints inside helices (except for the two
% N-terminal and C-terminal residues of the helix) according to Hovmoeller
% et al.

switch res
    case 'P'
        phi = 4*(rand-0.5) - 61.0;
        psi = 4*(rand-0.5) - 36.5;
    case 'G'
        phi = 4*(rand-0.5) - 59.1;
        psi = 4*(rand-0.5) - 42.4;
    otherwise
        phi = 4*(rand-0.5) - 63.8;
        psi = 4*(rand-0.5) - 41.1;
end
phi = pi*phi/180;
psi = pi*psi/180;

function allowed = test_Rama(poi,backbone,Rama_res,sequence)

[phi,psi,~] = dihedrals(poi,backbone);

% fprintf(1,'C-terminal anchor omega: %4.1f\n',180*omega/pi);
% fprintf(1,'C-terminal anchor phi  : %4.1f\n',180*phi/pi);
% fprintf(1,'C-terminal anchor psi  : %4.1f\n',180*psi/pi);

phipoi = 181 + round(180*phi/pi);
if phipoi < 1, phipoi = 1; end
if phipoi >360, phipoi = 360; end
psipoi = 181 + round(180*psi/pi);
if psipoi < 1, psipoi = 1; end
if psipoi >360, psipoi = 360; end

switch sequence(poi)
    case 'P'
        allowed = Rama_res.allowed_P(psipoi,phipoi);
    case 'G'
        allowed = Rama_res.allowed_G(psipoi,phipoi);
    otherwise
        allowed = Rama_res.allowed_gen(psipoi,phipoi);
end

function [backbone,success] = try_fix_rama(k0,backbone,Rama_res,sequence,direction)

step0 = 5;
step = direction*step0*pi/180; % step size in Ca-Ca rotation
k = k0;
backbone_f = backbone;

success = false;

while k > 1
    CA_back = backbone(4*k-6,:);
    C_back = backbone(4*k-5,:);
    O_back = backbone(4*k-4,:);
    N = backbone(4*k-3,:);
    CA = backbone(4*k-2,:);
    allowed = test_Rama(k,backbone,Rama_res,sequence);
    if ~allowed
        success = false;
        rotax = CA_back - CA;
        transmat = affine('rotn',step,rotax);
        for kr = 1:360/step0 
            C_back = affine_trafo_vector(C_back - CA,transmat) + CA;
            O_back = affine_trafo_vector(O_back - CA,transmat) + CA;
            N = affine_trafo_vector(N - CA,transmat) + CA;
            backbone_f(4*k-5,:) = C_back;
            backbone_f(4*k-4,:)  = O_back;
            backbone_f(4*k-3,:) = N;
            allowed = test_Rama(k,backbone_f,Rama_res,sequence);
            if allowed
                break
            end
        end
        if ~allowed
            return
        else
            success = true;
        end
        backbone(4*k-5,:) = C_back;
        backbone(4*k-4,:)  = O_back;
        backbone(4*k-3,:) = N;
        k = k - 1;
    else
        success = true;
        return
    end
end

function errcode = test_constraints(backbone,restrain,options,kend)

p_model = 1;
errcode = 0;

for k = 2:kend
    if ~isempty(restrain(k-1).label) 
        % make spin label coordinate
        x= backbone(4*k-3,:) - backbone(4*k-2,:); % x axis is along C_alpha-N bond
        x=x/norm(x);    % unit vector along x
        yp=backbone(4*k-1,:) - backbone(4*k-2,:); % y axis is in the plane spanned by x axis and C-Ca bond
        yp=yp/norm(yp);
        z=cross_rowvec(x,yp); % z axis is perpendicular on xy plane
        z=z/norm(z);
        y=cross_rowvec(z,x); % real (corrected) y axis 
        dircos=[x;y;z];
        Rp=dircos; % rotation matrix for conversion to standard frame
        restrain(k-1).xyz = restrain(k-1).label*Rp + backbone(4*k-2,:);
        p_beacon = 1;
        for kr = 1:length(restrain(k-1).r_beacon)
            r = norm(restrain(k-1).xyz-restrain(k-1).r_beacon(kr).xyz); 
            switch restrain(k-1).r_beacon(kr).type
                case 'Gaussian'
                    restrain(k-1).r_beacon(kr).p = ...
                        exp(-((r-restrain(k-1).r_beacon(kr).par1)/(sqrt(2)*restrain(k-1).r_beacon(kr).par2))^2);
                    p_beacon = p_beacon * restrain(k-1).r_beacon(kr).p;
                case 'bounds'
                    if r >= restrain(k-1).r_beacon(kr).par1 && ...
                            r <= restrain(k-1).r_beacon(kr).par2 % bounds violated
                       restrain(k-1).r_beacon(kr).p = 1;
                       p_beacon = 1;
                    end                        
                otherwise
                    error('MMMx:loop_model:unknownRestraintType','Restraint type %s not known',restrain(k-1).r_beacon(kr).type);
            end
        end
        p_model = p_beacon*p_model;
        if p_model < min_prob
            errcode = 5;
            return
        end
        p_intern = 1;
        for kr = 1:length(restrain(k-1).r_intern)
            site = restrain(k-1).r_intern(kr).site;
            r = norm(restrain(k-1).xyz-restrain(site).xyz); 
            switch restrain(k-1).r_intern(kr).type
                case 'Gaussian'
                    restrain(k-1).r_intern(kr).p = ...
                        exp(-((r-restrain(k-1).r_intern(kr).par1)/(sqrt(2)*restrain(k-1).r_intern(kr).par2))^2);
                    p_intern = p_intern * restrain(k-1).r_intern(kr).p;                 
                case 'bounds'
                    if r >= restrain(k-1).r_intern(kr).par1 && ...
                            r <= restrain(k-1).r_intern(kr).par2 % bounds violated
                        restrain(k-1).r_intern(kr).p = 1;
                        p_intern = 1;
                    end                        
                otherwise
                    error('MMMx:loop_model:unknownRestraintType','Restraint type %s not known',restrain(k-1).r_beacon(kr).type);
            end
        
        end
        p_model = p_intern*p_model;
        if p_model < min_prob
            errcode = 5;
            return
        end
        p_oligomer = 1;
        for kr = 1:length(restrain(k-1).oligomer)
            r = 2*sqrt(sum(restrain(k-1).xyz(1:2).^2))*sin(pi/restrain(k-1).oligomer(kr).n); 
            switch restrain(k-1).oligomer.type
                case 'Gaussian'
                    restrain(k-1).oligomer(kr).p = ...
                        exp(-((r-restrain(k-1).oligomer(kr).par1)/(sqrt(2)*restrain(k-1).oligomer(kr).par2))^2);
                    p_oligomer = p_oligomer * restrain(k-1).oligomer(kr).p;
                     
                case 'bounds'
                    if r >= restrain(k-1).oligomer(kr).par1 && ...
                            r <= restrain(k-1).oligomer.par2 % bounds violated
                        restrain(k-1).oligomer(kr).p = 1;
                        p_oligomer = 1;
                    end                        
                otherwise
                    error('MMMx:loop_model:unknownRestraintType','Restraint type %s not known',restrain(k-1).oligomer(kr).type);
            end
        end
        p_model = p_oligomer*p_model;
        if p_model < options.min_prob
            errcode = 5;
            return
        end
        p_depth = 1;
        for kr = 1:length(restrain(k-1).depth)
            switch restrain(k-1).depth(kr).site
                case 'CA'
                    z = abs(backbone(4*k-6,3));
                case 'label'
                    z = abs(restrain(k-1).xyz(3));
                otherwise
                    error('MMMx:loop_model:unknownRestraintModifier','Depth restraint site %s not known',restrain(k-1).depth(kr).site);
            end
            switch restrain(k-1).depth.type
                case 'Gaussian'
                    restrain(k-1).depth(kr).p = ...
                        exp(-((z-restrain(k-1).depth(kr).par1)/(sqrt(2)*restrain(k-1).depth(kr).par2))^2);
                    p_depth = p_depth * restrain(k-1).depth(kr).p;                 
                case 'bounds'
                    if z >= restrain(k-1).depth(kr).par1 && ...
                            z <= restrain(k-1).depth.par2 % bounds violated
                        restrain(k-1).depth(kr).p = 1;
                        p_depth = 1;
                    else
                        p_depth = 0;
                    end                       
                otherwise
                    error('MMMx:loop_model:unknownRestraintType','Restraint type %s not known',restrain(k-1).depth(kr).type);
            end
        end
        p_model = p_depth*p_model;
        if p_model < options.min_prob
            errcode = 5;
            return
        end
    end
end