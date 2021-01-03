function [coor,errcode,restrain,p_model,k] = loop_model_reverse(sequence, anchorC, anchorCn, prot_xyz, restrain, Rama_res, rescodes, n_res, options)
% [coor,errcode,restrain,p_model,k] = loop_model_reverse(sequence, anchorC, anchorCn, prot_xyz, restrain, Rama_res, rescodes, n_res, options)
% Generates a loop model anchored at the N terminus of an existing protein 
% or peptide structure, conforming to residue-specific Ramachandran plots
% backbone based on: H. Sugeta, T. Miyazawa, Biopolymers, 1967, 5, 673-679.
% residue-specific Ramachandran plots from: S. Hovmöller, T. Zhou, T.
% Ohlson, Acta Cryst. D, 2002, 58, 768-776.
%
% Input:
% sequence  peptide sequence in single-letter code, must contain C-terminal
%           anchor
% anchorC   backbone coordinates of C-terminal anchor residue in the order
%           N, CA, C, O (Cartesian in Å) 
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
% rescodes  MMM-internal residue codes for the sequence, must contain
%           C-terminal anchor
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
%           .max_res            maximum attempts for a single-residue step,
%                               defaults to 100
%
% Output:
% coor      Cartesian backbone coordinates of the loop in the order N, CA,
%           C, O; empty if not successful
% errcode   error code for unsuccessful attempts, 0 for success
%           2   loop had a self clash
%           4   loop clashed with protein
%           5   restraint violation
% restrain  restraint variable with diagnostic information
%           .xyz    simulated label position
%           .r_beacon(l).p  fulfillment probability
%           .r_intern(l).p  fulfillment probability
%           .oligomer.p  fulfillment probability
%           .depth.p  fulfillment probability
% p_model   cumulative model probability from restraints
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
beta_psi_LB = 128;
beta_psi_UB = 145;

% PPII helix region
PPII_phi = -75;
PPII_psi = 160;

% cis peptide frequency according to Jena library of Biological Macromolecules

cis_Xaa_Pro = 0.05;
cis_Xaa_non_Pro = 0.0003;

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

ngap = length(sequence)-1;

backbone = zeros(4*(ngap+1),3); 
backbone(end-3:end,:) = anchorC; 

phivec = zeros(1,ngap);
psivec = zeros(1,ngap);

% bootstrapping, phi of C terminal anchor not yet defined
N_forth = anchorCn(1,:);
N = anchorC(1,:);
CA = anchorC(2,:);
C = anchorC(3,:);
psi = dihedral_fast(N,CA,C,N_forth);

poi = round(rand*Rama_res.me{rescodes(end)}+0.5);
ephi=Rama_res.ephi{rescodes(end)}(poi);
epsi=Rama_res.epsi{rescodes(end)}(poi);

[~,poi] = min(abs(epsi-psi)); % find close match of C-terminal anchor psi
phi = ephi(poi); % set a corresponding psi value

% local frame of Calpha of the N-terminal anchor residue
x = anchorC(1,:)-anchorC(2,:);
bondx = norm(x);
x = x/norm(x);
yp = anchorC(2,:)-anchorC(3,:);
yp = yp/norm(yp);
z = cross_rowvec(x,yp);
z = z/norm(z);
y = cross_rowvec(z,x);
y = y/norm(y);
A = [x;y;z];
% transformation matrix into that local frame
A = A';
% coordinates of Calpha of the N-terminal anchor residue
acoor = anchorC(2,:);
acoor = acoor';
acoor = acoor + A*[bondx;0;0]; % coordinates of N of the C-terminal anchor residue

ctau=cos(-phi); % sign changed for 'backward' extension
stau=sin(-phi); % sign changed for 'backward' extension
A23=[-cfi3,-sfi3,0;sfi3*ctau,-cfi3*ctau,-stau;sfi3*stau,-cfi3*stau,ctau];
A=A*A23;
acoor=acoor+A*b3vec'; % coordinates of C of first missing residue
% decide whether cis or trans peptide
so = sot;
co = cot;
dice = rand;
if upper(sequence(1))~='P'
    if dice < cis_Xaa_non_Pro
        so = soc; co = coc;
    end
else
    if dice < cis_Xaa_Pro
        so = soc; co = coc;
    end
end
A31=[-cfi2,-sfi2,0;sfi2*co,-cfi2*co,-so;sfi2*so,-cfi2*so,co];
A=A*A31; % local frame at C of first missing residue

reboot.k = length(sequence);
reboot.counter = reboots;
reboot.A = A;
reboot.acoor = acoor;
reboot.p_model = p_model;
reboot.tested = 0;

reboot_thresh = options.min_prob^(1/n_res);
tested_restraints = 0;

% loop generation
k = length(sequence)-1;
while k >= 1
    updated = false;
    % resnum = residues(k);
    if restrain(k).secondary == 3
        [rphi,rpsi] = get_phi_psi_in_helix(sequence(k));
    elseif restrain(k).secondary == 4
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
    if upper(sequence(k))~='P'
        if dice < cis_Xaa_non_Pro
            so = soc; co = coc;
        end
    else
        if dice < cis_Xaa_Pro
            so = soc; co = coc;
        end
    end
    if restrain(k).cis
        so = soc;
        co = coc;
    end
    if restrain(k).secondary
        att_sec = 1;
        if restrain(k).secondary == 1 % alpha-helix
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
        elseif restrain(k).secondary == 2 % beta-sheet
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
    backbone(4*k-1,:)=acoor'; % coordinates of C
    acoor=acoor+A*b2vec'; % coordinates of C_alpha
    backbone(4*k-2,:)=acoor';
    ctau=cos(-rpsi); % sign changed for 'backward' prediction
    stau=sin(-rpsi); % sign changed for 'backward' prediction
    A12=[-cfi1,-sfi1,0;sfi1*ctau,-cfi1*ctau,-stau;sfi1*stau,-cfi1*stau,ctau];
    A=A*A12;
    acoor=acoor+A*b1vec'; % coordinates of N
    backbone(4*k-3,:)=acoor';
    if ~isempty(restrain(k).label)
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
        restrain(k).xyz = restrain(k).label*Rp + backbone(4*k-2,:);
        p_beacon = 1;
        for kr = 1:length(restrain(k).r_beacon)
            r = norm(restrain(k).xyz-restrain(k).r_beacon(kr).xyz); 
            switch restrain(k).r_beacon(kr).type
                case 'Gaussian'
                    restrain(k).r_beacon(kr).p = ...
                        exp(-((r-restrain(k).r_beacon(kr).par1)/(sqrt(2)*restrain(k).r_beacon(kr).par2))^2);
                    p_beacon = p_beacon * restrain(k).r_beacon(kr).p;
                    tested_restraints = tested_restraints + 1;
                    updated = true;
                case 'bounds'
                    if r < restrain(k).r_beacon(kr).par1 || ...
                            r > restrain(k).r_beacon(kr).par2 % bounds violated
                        restrain(k).r_beacon(kr).p = 0;
                        p_beacon = 0;
                    end                        
                otherwise
                    error('MMMx:loop_model_reverse:unknownRestraintType','Restraint type %s not known',restrain(k).r_beacon(kr).type);
            end
        end
        p_model = p_beacon*p_model;
        if p_model < options.min_prob
            if reboot.counter > 0
                reboot.counter = reboot.counter - 1;
                A = reboot.A;
                acoor = reboot.acoor;
                k = reboot.k - 1;
                tested_restraints = reboot.tested;
                p_model = reboot.p_model;
                continue
            end
            errcode = 5;
            coor = [];
%                 fprintf(2,'Model rejected at residue %i by too low restraint fulfillment (beacon).\n',k);
            return
        end
        p_intern = 1;
        for kr = 1:length(restrain(k).r_intern)
            site = restrain(k).r_intern(kr).site;
            r = norm(restrain(k).xyz-restrain(site).xyz); 
            switch restrain(k).r_intern(kr).type
                case 'Gaussian'
                    restrain(k).r_intern(kr).p = ...
                        exp(-((r-restrain(k).r_intern(kr).par1)/(sqrt(2)*restrain(k).r_intern(kr).par2))^2);
                    p_intern = p_intern * restrain(k).r_intern(kr).p;     
                    tested_restraints = tested_restraints + 1;
                    updated = true;
                case 'bounds'
                    if r < restrain(k).r_intern(kr).par1 || ...
                            r > restrain(k).r_intern(kr).par2 % bounds violated
                       restrain(k).r_intern(kr).p = 0;
                       p_intern = 0;
                    end                        
                otherwise
                    error('MMMx:loop_model_reverse:unknownRestraintType','Restraint type %s not known',restrain(k).r_beacon(kr).type);
            end
        end
        p_model = p_intern*p_model;
        if p_model < options.min_prob
            if reboot.counter > 0
                reboot.counter = reboot.counter - 1;
                A = reboot.A;
                acoor = reboot.acoor;
                k = reboot.k - 1;
                tested_restraints = reboot.tested;
                p_model = reboot.p_model;
                continue
            end
            errcode = 5;
            coor = [];
            return
        end
        p_oligomer = 1;
        for kr = 1:length(restrain(k).oligomer)
            r = 2*sqrt(sum(restrain(k).xyz(1:2).^2))*sin(pi/restrain(k).oligomer(kr).n); 
            switch restrain(k).oligomer.type
                case 'Gaussian'
                    restrain(k).oligomer(kr).p = ...
                        exp(-((r-restrain(k).oligomer(kr).par1)/(sqrt(2)*restrain(k).oligomer(kr).par2))^2);
                    p_oligomer = p_oligomer * restrain(k).oligomer(kr).p;
                    tested_restraints = tested_restraints + 1; 
                    updated = true;
                case 'bounds'
                    if r < restrain(k).oligomer(kr).par1 || ...
                            r > restrain(k).oligomer.par2 % bounds violated
                        restrain(k).oligomer(kr).p = 0;
                        p_oligomer = 0;
                    end                       
                otherwise
                    error('MMMx:mk_loop_model_reverse:unknownRestraintType','Restraint type %s not known',restrain(k).oligomer(kr).type);
            end
        end
        p_model = p_oligomer*p_model;
        if p_model < options.min_prob
            if reboot.counter > 0
                reboot.counter = reboot.counter - 1;
                A = reboot.A;
                acoor = reboot.acoor;
                k = reboot.k - 1;
                tested_restraints = reboot.tested;
                p_model = reboot.p_model;
                continue
            end
            errcode = 5;
            coor = [];
            return
        end
        p_depth = 1;
        for kr = 1:length(restrain(k).depth)
            switch restrain(k).depth(kr).site
                case 'CA'
                    z = abs(backbone(4*k-2,3));
                case 'label'
                    z = abs(restrain(k).xyz(3));
                otherwise
                    error('MMMx:loop_model_reverse:unknownRestraintModifier','Depth restraint site %s not known',restrain(k).depth(kr).site);
            end
            switch restrain(k).depth.type
                case 'Gaussian'
                    restrain(k).depth(kr).p = ...
                        exp(-((z-restrain(k).depth(kr).par1)/(sqrt(2)*restrain(k).depth(kr).par2))^2);
                    p_depth = p_depth * restrain(k).depth(kr).p;  
                    tested_restraints = tested_restraints + 1;
                    updated = true;
                case 'bounds'
                    if z >= restrain(k).depth(kr).par1 && ...
                            z <= restrain(k).depth.par2 % bounds not violated
                        restrain(k).depth(kr).p = 1;
                        p_depth = 1;
                    else
                        p_depth = 0;
                    end                        
                otherwise
                    error('MMMx:loop_model_reverse:unknownRestraintType','Restraint type %s not known',restrain(k).depth(kr).type);
            end
        end
        p_model = p_depth*p_model;
        if p_model < options.min_prob
            if reboot.counter > 0
                reboot.counter = reboot.counter - 1;
                A = reboot.A;
                acoor = reboot.acoor;
                k = reboot.k - 1;
                tested_restraints = reboot.tested;
                p_model = reboot.p_model;
                continue
            end
            errcode = 5;
            coor = [];
            return
        end
    end
    
    ctau=cos(-rphi); % sign changed for 'backward' prediction
    stau=sin(-rphi); % sign changed for 'backward' prediction
    A23=[-cfi3,-sfi3,0;sfi3*ctau,-cfi3*ctau,-stau;sfi3*stau,-cfi3*stau,ctau];
    A=A*A23;
    acoor=acoor+A*b3vec'; % coordinates of next C
    A31=[-cfi2,-sfi2,0;sfi2*co,-cfi2*co,-so;sfi2*so,-cfi2*so,co];
    A=A*A31;    
    if p_model^(1/tested_restraints) > reboot_thresh && updated % set reboot point, if model is 'better than expected' at this poiunt
        reboot.k = k;
        reboot.counter = reboots;
        reboot.A = A;
        reboot.acoor = acoor;
        reboot.p_model = p_model;
        reboot.tested = tested_restraints;
    end
    k = k - 1;
end

% if the model survived up to this point, it has sufficient probability of
% fulfilling all constraints (or flag 'constrained' was set to false)

% generate backbone oxygen coordinates
for k = 1:length(sequence)-1
    backbone = add_O(k,backbone);
end

% check for clashes with protein
pair_dist = get_all_pair_dist(backbone(1:4*(length(sequence)-2),:),prot_xyz);
min_dist = min(min(pair_dist));
if min_dist < options.clash_threshold_lp
    errcode = 4;
    p_model = 0;
    coor = [];
    return
end

% self-clash test
k = 1;
clash = false;
while ~clash && k <= length(sequence)-2
    min_dist = get_min_pair_dist(k,backbone);
    if min_dist < options.clash_threshold
        clash = true;
    end
    k = k + 1;
end

if clash
    errcode = 2;
    coor = [];
    p_model = 0;
    return
end

coor = backbone;

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