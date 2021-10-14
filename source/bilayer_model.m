function [entity,error] = bilayer_model(entity,conformer,type,oriented,logfid)
%
% BILAYER_MODEL Compute lipid bilayer plane and thickness that minimizes
%               free energy of protein immersion 
%
%   [transform,energy] = BILAYER_MODEL(entity,conformer)
%   tansformation parameters and free energy of protein immersion for a
%   bilayer model
%
%   based on: L. Adamian, V. Nanda, W. F. De Grado, J. Liang,
%             Proteins, 2005, 59:496-509
%
% INPUT
% entity        protein structure in MMMx:atomic format
% conformer     optional conformer number, defaults to 1
% type          membrane protein type, 'bundle' or 'barrel'
% oriented      flag indicating that the z axis of the input is the bilayer
%               normal, defaults to false (orientation is fitted)
% logfid        optional file identifier of a logfile, if missing, output
%               is directed to the Matlab console, use empty logfid to
%               suppress output
%
% OUTPUT
% bilayer       bilayer definition, struct with fields
%               .theta      polar angle theta of membrane normal
%               .phi        polar angle phi of membrane normal
%               .z          shift of membrane center plane along membrane
%                           normal
%               .thickness  bilayer thickness 
% energy        free energy of immersion [kcal/mol]
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2021: Gunnar Jeschke

if ~exist('conformer','var') || isempty(conformer)
    conformer = 1;
end

if ~exist('type','var') || isempty(type)
    type = 'bundle';
end

if ~exist('oriented','var') || isempty(oriented)
    oriented = false;
end

if ~exist('logfid','var')
    logfid = 1;
elseif isempty(logfid)
    logfid = 3;
end

maxit = 10; % maximum number of iterations

energy_convergence = 1; % 1 kcal/mol energy convergence criterion
z_convergence = 0.2; % 0.2 Angstroem z shift convergence criterion
theta_convergence = 2; % 2 degeree theta convergence criterion
phi_convergence = 5; % 5 degree phi convergence criterion

% default parameters

probe_radius = 1.9; % based on L. Adamian, V. Nanda, W. F. De Grado, J. Liang
                  % Proteins, 2005, 59:496-509
bilayer.width = 60; % initial width of membrane slab
bilayer.thickness = 35; % initial layer thickness
bilayer.z = 0; % initial center shift
bilayer.theta = 0; % initial membrane normal angle theta
bilayer.phi = 0; % initial membrane normal angle phi

[TM_addresses,propensities] = get_TM_residues(entity);

switch type
    case 'bundle'
        propensities = propensities(:,1:2);
    case 'barrel'
        propensities = propensities(:,3:4);
    otherwise
        error = sprintf('Unknown type ''%s'' of membrane protein',type); 
        return
end

[access,error] = accessibilities(entity,conformer,TM_addresses,probe_radius);

if ~isempty(error)
    return
end

accessibility = access(:,2);

energy0 = immersion_energy(entity,conformer,TM_addresses,propensities,accessibility,bilayer);

fprintf(logfid,'\nBilayer fit. Initial immersion free energy is %6.1f kcal/mol\n',energy0);

it = 0;
converged = false;
bilayer0 = bilayer;
while it < maxit && ~converged
    it = it + 1;
    [bilayer,energy] = fit_bilayer(entity,conformer,TM_addresses,propensities,accessibility,bilayer0,oriented);
    
    fprintf(logfid,'Iteration %i. Immersion free energy is %6.1f kcal/mol\n',it,energy);
    fprintf(logfid,'thickness   : %4.1f %s\n',bilayer.thickness,char(197));
    fprintf(logfid,'center shift: %4.1f %s\n',bilayer.z,char(197));
    fprintf(logfid,'altitude    : %4.1f %s\n',180*bilayer.theta/pi,char(176));
    fprintf(logfid,'azimuth     : %4.1f %s\n',180*bilayer.phi/pi,char(176));
    
    delta_z = abs(bilayer0.z - bilayer.z);
    delta_theta = 180*abs(bilayer0.theta - bilayer.theta);
    delta_phi = 180*abs(bilayer0.phi - bilayer.phi);

    if abs(energy - energy0) < energy_convergence
        converged = true;
    end
    if delta_z < z_convergence...
            && delta_theta < theta_convergence...
            && delta_phi < phi_convergence
        converged = true;
    end
    energy0 = energy;
    transmat1 = affine('Euler',[bilayer.phi,bilayer.theta,0]);
    transmat2 = affine('translation',[0,0,-bilayer.z]);
    entity.xyz = affine_coor_set(entity.xyz,{transmat1,transmat2});

end
if converged
    fprintf(logfid,'Bilayer fit converged\n');
else
    fprintf(logfid,'Bilayer fit stopped after %i iterations\n',maxit);
end

function [TM_addresses,propensities] = get_TM_residues(entity)

defs = load('monomer_attributes.mat');
aa = defs.monomers.aa;
aa_tags = defs.monomers.aa_tags;

TM_poi = 0;
propensities = zeros(50000,2);
TM_addresses = cell(10000,1);

chains = fieldnames(entity);

for kc = 1:length(chains)
    chain = chains{kc};
    if isstrprop(chain(1),'upper') % chain fields start with a capital
        residues = fieldnames(entity.(chain));
        for kr = 1:length(residues) % expand over all residues
            residue = residues{kr};
            if strcmp(residue(1),'R') % these are residue fields
                aa_index = tag2id(entity.(chain).(residue).name,upper(aa_tags));
                if ~isempty(aa_index)
                    TM_poi = TM_poi + 1;
                    propensities(TM_poi,1) = aa.TMLIPH(aa_index);
                    propensities(TM_poi,2) = aa.TMLIPC(aa_index);
                    propensities(TM_poi,3) = aa.TMBH(aa_index);
                    propensities(TM_poi,4) = aa.TMBC(aa_index);
                    TM_addresses{TM_poi} = sprintf('%s.%s',chain,residue);
                end
            end
        end
    end
end

propensities = propensities(1:TM_poi,:);
TM_addresses = TM_addresses(1:TM_poi);

function [access,error] = accessibilities(entity,conformer,TM_addresses,probe_radius)

defs = load('monomer_attributes.mat');
aa = defs.monomers.aa;
aa_tags = defs.monomers.aa_tags;

error = '';

tolerance = 1.2; % tolerance for surface distance in ?

chemistry = load('element_attributes.mat');

dospath = which_third_party_module('msms');
if isempty(dospath)
    access = [];
    error = 'Bilayer computation requires MSMS by M. Sanner';
    return
end

elements = entity.elements;
coor = entity.xyz;
vdW=zeros(size(elements));
for k = 1:length(elements)
    if elements(k)>0 && elements(k) < length(chemistry.pse.vdW)
        vdW(k) = chemistry.pse.vdW(elements(k));
    end
end

density = round(5000/length(elements));
if density>3, density = 3; end
if density<1, density = 1; end
dstring = sprintf(' -density %i',density);


radstr = sprintf(' -probe_radius %4.1f',probe_radius);

outfile = 'tmp.xyzr';
outfile = fullfile(pwd,outfile);
fid = fopen(outfile,'w');
if fid == -1
    access = [];
    return;
end
for k = 1:length(vdW)
    if vdW(k)>0 
        fprintf(fid,'%12.3f%12.3f%12.3f%6.2f',coor(k,1),coor(k,2),coor(k,3),vdW(k));
        if k < length(vdW)
            fprintf(fid,'\n');
        end
    end
end
fclose(fid);

msmsfile = 'msms_tmp';
msmsfile = fullfile(pwd,msmsfile);
cmd=[dospath ' -if ' outfile ' -of ' msmsfile dstring radstr];
[s,~] = dos(cmd);
if s ~= 0
    error = 'MSMS computation of solvent-accessible surface failed';
    access = [];
    return
end

[tri,x,y,z] = rd_msms(msmsfile);

delete([msmsfile '.face']);
delete([msmsfile '.vert']);


[mt,~]=size(tri);

midpoints = zeros(mt,3);
area = zeros(1,mt);
for k=1:mt
    vert = [x(tri(k,:))',y(tri(k,:))',z(tri(k,:))'];
    a = norm(vert(1,:)-vert(2,:));
    b = norm(vert(1,:)-vert(3,:));
    c = norm(vert(2,:)-vert(3,:));
    midpoints(k,:) = mean(vert);
    s = (a+b+c)/2;
    area(k) = sqrt(s*(s-a)*(s-b)*(s-c));
end


m = length(TM_addresses);

access=zeros(10000,2);
poi=0;
for k = 1:m
    adr_parts = split(TM_addresses{k},'.');
    [xyz,elements] = get_residue_atoms(entity,adr_parts{1},adr_parts{2},conformer);
    tag = entity.(adr_parts{1}).(adr_parts{2}).name;
    aa_index = tag2id(tag,upper(aa_tags));
    if ~isempty(aa_index)
        A = aa.A(aa_index);
    else
        A = 0;
    end
    mask = zeros(1,mt);
    for kkk = 1:length(elements)
        if elements(kkk) > 0 && elements(kkk) < length(chemistry.pse.vdW)
            vdW = chemistry.pse.vdW(elements(kkk));
            rad = (vdW+tolerance)^2;
            atom = repmat(xyz(kkk,:),mt,1);
            dvec = midpoints - atom;
            dist = sum(dvec.^2,2) - rad;
            dist = dist<=0;
            mask = mask + dist';
        end
    end
    myarea = sum(area(mask > 0));
    if A ~= 0
        relacc = myarea/A;
    else
        relacc = -1;
    end
    poi = poi+1;
    access(poi,1) = myarea;
    access(poi,2) = relacc;
end
access=access(1:poi,:);


function [tri,x,y,z,info] = rd_msms(filename)
% function [tri,x,y,z] = rd_msms(filename)
% 
% Reads the triangle definitions tri and vertices (x,y,z) that define the
% solvent excluded surface as computed with the MSMS program 
% by Michel F. Sanner 
% (see: http://www.scripps.edu/~sanner/html/msms_home.html)
%
% the proper citation is:
% Sanner, M. F., Spehner, J.-C., Olson, A. J. (1996) Reduced surface: an
% efficient way to compute molecular surfaces. Biopolymers, 38(3), 305-320.
%
% the program MSMS is not (and must not be) bundled with MMM, users have to
% obtain it for themselves or for their institution
%
% files filname.face and filename.vert must exist, otherwise the
% corresponding variables are empty, info provides parameters
%
% G. Jeschke, 2009

tri=[];

info.header='';
info.faces=[];
info.spheres=[];
info.density=[];
info.probe_r=[];

fid=fopen([filename '.face'],'r');

if fid~=-1
    tline = fgetl(fid);
    info.header = strtrim(tline(2:end));
    fgetl(fid);
    info.faces = fscanf(fid,'%i',1);
    info.spheres = fscanf(fid,'%i',1);
    info.density = fscanf(fid,'%f',1);
    info.probe_r = fscanf(fid,'%f',1);
    fgetl(fid);
    tri=zeros(info.faces,3);
    for k = 1:info.faces
        tline = fgetl(fid);
        tri(k,:) = sscanf(tline,'%i',3);
    end
end

fclose(fid);

fid=fopen([filename '.vert'],'r');

if fid ~= -1
    fgetl(fid);
    fgetl(fid);
    info.vertices = fscanf(fid,'%i',1);
    coor = zeros(info.vertices,3);
    fgetl(fid);
    for k = 1:info.vertices
        tline = fgetl(fid);
        coor(k,:) = sscanf(tline,'%f',3);
    end
end

fclose(fid);

x=coor(:,1)';
y=coor(:,2)';
z=coor(:,3)';

function [xyz,elements] = get_residue_atoms(entity,chain,residue,conformer)

xyz = zeros(100,3);
elements = zeros(100,1);
atnum = 0;

atoms = fieldnames(entity.(chain).(residue));
for ka = 1:length(atoms) % expand over all atoms
    atom = atoms{ka};
    if isstrprop(atom(1),'upper') % these are atom fields
        index = entity.(chain).(residue).(atom).tab_indices(conformer);
        atnum = atnum + 1;
        xyz(atnum,:) = entity.xyz(index,:);
        elements(atnum) = entity.elements(index);
    end
end
xyz = xyz(1:atnum,:);
elements = elements(1:atnum);

function energy = immersion_energy(entity,conformer,TM_addresses,propensities,accessibility,bilayer)
%
% hydrophobic energy contribution of residues an insertion into a bilayer
% the bilayer normal is assumed to be the z axis
% based on: L. Adamian, V. Nanda, W. F. De Grado, J. Liang,
%           Proteins, 2005, 59:496-509
%

transmat = affine('Euler',[bilayer.phi,bilayer.theta,0]);

RT = 8.314*(273.15+27)/4.1868; % RT at 27 ?C in kcal/mol

m = length(TM_addresses);

succ=0;
energy = 0;
for k = 1:m
    adr_parts = split(TM_addresses{k},'.');
    chain = adr_parts{1};
    residue = adr_parts{2};
    if isfield(entity.(chain).(residue),'CA')
        index = entity.(chain).(residue).CA.tab_indices(conformer);
        coor = [entity.xyz(index,:) 1];
        coor = transmat*coor';
        coor(3) = coor(3) - bilayer.z;
        position = abs(coor(3))/bilayer.thickness;
        propensity = 0;
        if position <= 0.5 % headgroup layer
            propensity = propensities(k,1);
        end
        if position < 0.25 % core
            propensity = propensities(k,2);
        end
        energy = energy + accessibility(k)*propensity;
        succ=1;
    end
end

if ~succ 
    energy = []; 
else
    energy = RT*energy;
end

function [bilayer,energy] = fit_bilayer(entity,conformer,TM_addresses,propensities,accessibility,bilayer,oriented)

z_grid = linspace(-10,10,21);
d_grid = linspace(25,50,21);
min_energy = 1e10;

opt_z = bilayer.z;
opt_d = bilayer.thickness;
c_bilayer = bilayer;

options = optimset('fminsearch');
options = optimset(options,'TolFun',0.05,'TolX',0.05);

for k = 1:length(z_grid)
    c_bilayer.z = z_grid(k);
    for kk=1:length(d_grid)
        c_bilayer.thickness = d_grid(kk);
        energy = immersion_energy(entity,conformer,TM_addresses,propensities,accessibility,c_bilayer);
        if energy < min_energy
            min_energy = energy;
            opt_z = c_bilayer.z;
            opt_d = c_bilayer.thickness;
        end
    end
end

c_bilayer.z = opt_z;
c_bilayer.thickness = opt_d;
c_bilayer.theta = 0;
c_bilayer.phi = 0;

if oriented
    v0 = [opt_z opt_d];
    [v,energy] = fminsearch(@energy_fit_oriented,v0,options,entity,conformer,TM_addresses,propensities,accessibility,c_bilayer);
    bilayer.z = v(1);
    bilayer.thickness = v(2);
else
    th_grid = linspace(-pi/4,pi/4,21);
    phi_grid = linspace(-pi,pi,41);
    
    min_energy = 1e10;
    opt_theta = bilayer.theta;
    opt_phi = bilayer.phi;
    
    for k = 1:length(th_grid)
        c_bilayer.theta = th_grid(k);
        for kk = 1:length(phi_grid)
            c_bilayer.phi = phi_grid(kk);
            energy = immersion_energy(entity,conformer,TM_addresses,propensities,accessibility,c_bilayer);
            if energy < min_energy
                min_energy = energy;
                opt_theta = c_bilayer.theta;
                opt_phi = c_bilayer.phi;
            end
        end
    end
    
    c_bilayer.theta = opt_theta;
    c_bilayer.phi = opt_phi;
    v0 = [opt_z opt_d opt_theta opt_phi];
    [v,energy] = fminsearch(@energy_fit,v0,options,entity,conformer,TM_addresses,propensities,accessibility,c_bilayer);
    bilayer.z = v(1);
    bilayer.thickness = v(2);
    bilayer.theta = v(3);
    bilayer.phi = v(4);
end




function energy = energy_fit(v,entity,conformer,TM_addresses,propensities,accessibility,bilayer)

bilayer.z = v(1);
bilayer.thickness = v(2);
bilayer.theta = v(3);
bilayer.phi = v(4);

energy = immersion_energy(entity,conformer,TM_addresses,propensities,accessibility,bilayer);

function energy = energy_fit_oriented(v,entity,conformer,TM_addresses,propensities,accessibility,bilayer)

bilayer.z = v(1);
bilayer.thickness = v(2);

energy = immersion_energy(entity,conformer,TM_addresses,propensities,accessibility,bilayer);