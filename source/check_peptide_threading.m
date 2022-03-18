function threaded = check_peptide_threading(coor_prot,coor_RNA,atomtags)
%
% CHECK_PEPTIDE_THREADING   Checks if a peptide is threaded through a hole
%                           in the RNA
%
%   threaded = CHECK_PEPTIDE_THREADING(coor_prot,coor_RNA,atomtags)
%   Returns a Boolean flag that indicates whether a peptide is close to
%   RNA atoms in all directions, a single residue is always
%   considered as not threaded
%
% INPUT
% coor_prot    [np,3] array of Cartesian protein coordinates
% coor_RNA     [nr,3] array of Cartesian RNA coordinates
% atomtags     atom tags corresponding to peptide coordinates, only CA and
%              O coordinates are used
%
% OUTPUT
% threaded     true, if the peptide is too close to RNA atoms in all
%              directions, false if not; empty if there is a mismatch
%              between the numbers of CA and O atoms
%

% This file is a part of MMMx. License is MIT (see LICENSE.md). 
% Copyright(c) 2022: Gunnar Jeschke

threshold = 3.6; % 3.6, threshold for close approach, optimized on PTBP1/EMCV-IRES
test_radius1 = 2.0; % 2.0, smaller test radius, optimized on PTBP1/EMCV-IRES
test_radius2 = 5.0; % 5.0, intermediate test radius, optimized on PTBP1/EMCV-IRES
test_radius3 = 6.0; % 6.0, larger test radius, optimized on PTBP1/EMCV-IRES

n0 = length(atomtags);
CA_coor = zeros(n0,3);
O_coor = zeros(n0,3);
nCA = 0;
nO = 0;
for k = 1:n0
    if strcmp(atomtags{k},'CA')
        nCA = nCA + 1;
        CA_coor(nCA,:) = coor_prot(k,:);
    end
    if strcmp(atomtags{k},'O')
        nO = nO + 1;
        O_coor(nO,:) = coor_prot(k,:);
    end
end
if nCA ~= nO 
    threaded = [];
    return
else
    CA_coor = CA_coor(1:nCA,:);
    O_coor = O_coor(1:nO,:);
end
if nCA < 2
    threaded = false;
    return
end

spr = 6;
[backbone,rung,normal] = mk_ribbon(spr,CA_coor,O_coor);

angi = 0:pi/50:2*pi-pi/50;
sth = sin(angi);
cth = cos(angi);
[nr,~] = size(coor_RNA);
approaches1 = zeros(1,length(backbone));
approaches2 = zeros(1,length(backbone));
approaches3 = zeros(1,length(backbone));
for k = 1:length(backbone)
    approach1 = zeros(1,length(angi));
    approach2 = zeros(1,length(angi));
    approach3 = zeros(1,length(angi));
    for a = 1:length(angi)
        spike = cth(a)*rung(k,:) + sth(a)*normal(k,:);
        point = backbone(k,:) + test_radius1*spike;
        diff = sqrt(sum((coor_RNA - repmat(point,nr,1)).^2,2));
        approach1(a) = min(diff);
        point = backbone(k,:) + test_radius2*spike;
        diff = sqrt(sum((coor_RNA - repmat(point,nr,1)).^2,2));
        approach2(a) = min(diff);
        point = backbone(k,:) + test_radius3*spike;
        diff = sqrt(sum((coor_RNA - repmat(point,nr,1)).^2,2));
        approach3(a) = min(diff);
    end
    approaches1(k) = max(approach1);
    approaches2(k) = max(approach2);
    approaches3(k) = max(approach3);
end
approaches = max([approaches1;min([approaches2;approaches3])]);
% Diagnostic plotting
figure(1); clf; hold on;
plot(approaches1,'b.');
plot(approaches2,'.','Color',[0,0.6,0]);
plot(approaches3,'r.');
plot(approaches,'k');
drawnow
threaded = min(approaches) < threshold;

function [backbone,rung,normal]=mk_ribbon(spr,Ca_coor,O_coor,sheetflag)
% function [backbone,rung,normal]=mk_ribbon(Ca_coor,O_coor,sheetflag)
%
% Create coordinates of a ribbon guide curve for a protein backbone segment,
% using a slightly modified version of the Carson-Bugg algorithm (J. Mol. 
% Graph. 1986, 4, 121-122), unit vectors in and normal to the ribbon face
% are also supplied
%
% uses: function bezierInterp
%
% spr                   segments per residue
% Ca_coor               coordinates of Calpha atoms
% O_coor                optional, coordinates of peptide oxygen atoms
% sheetflag             0 for loops and alpha-helices, 1 for sheets,
%                       defaults to 0
%
% backbone              coordinates of backbone guide points
% rung                  unit vectors within peptide plane, perpendicular to
%                       backbone, is empty, if O_coor is missing
% normal                unit vectors perpendicular to both peptide plane
%                       and backbone, is empty, if O_coor is missing
% 
% (c) G. Jeschke, 2009

backbone=[];
rung=[];
normal=[];

Tension=0; % assume cardinal splines with zero tension

if nargin<4
    sheetflag=0;
end

[m,n]=size(Ca_coor);

if m<2
    return;
end

guide_length=(m-1)*spr+1;
backbone=zeros(guide_length,3);
backbone0=[Ca_coor(1,:); Ca_coor; Ca_coor(m,:)];

if nargin>2
    rung0=zeros(m,n);
    dirsign=1;

    % run over all Calpha atoms
    for k=1:m-1
        a=Ca_coor(k+1,:)-Ca_coor(k,:);
        a=a/norm(a);
        b=O_coor(k,:)-Ca_coor(k,:);
        c=cross_rowvec(a,b);
        c=c/norm(c); % normal to peptide plane
        d=cross_rowvec(c,a);
        d=d/norm(d); % y axis in peptide plane
        % if sheetflag, dirsign=dirsign*-1; end; % flip side, if necessary
        if k>1
            dirsign=dirsign*sum(d.*d0);
        end
        rung0(k,:)=dirsign*d;
        d0=d;
    end
    rung0(m,:)=dirsign*d;
    normal=backbone;
    rung=backbone;
    rung0=[rung0(1,:); rung0; rung0(m,:)];    
end

% Calculate backbone spline interpolant
for k=1:m-1
    basnum=(k-1)*spr;
    backbone(basnum+1:basnum+spr+1,:)=cardinal_spline(backbone0(k:k+3,:),Tension,spr);
    if nargin>2
       rung(basnum+1:basnum+spr+1,:)=cardinal_spline(rung0(k:k+3,:),Tension,spr);    
    end
end

if nargin>2
    normal=zeros(size(rung));
    for k=1:guide_length-1
        e=backbone(k+1,:)-backbone(k,:);
        f=rung(k,:);
        g=cross_rowvec(e,f);
        normal(k,:)=g/norm(g);
        rung(k,:)=f/norm(f); %#ok<AGROW>
    end
    normal(guide_length,:)=normal(guide_length-1,:);
    rung(guide_length,:)=rung(guide_length,:)/norm(rung(guide_length,:));
end

if sheetflag
    backbone=smooth_backbone(backbone,spr);
end

function backbone=smooth_backbone(backbone,spr)
% Polynomial fitting (and smoothing) of the backbone curve

[m,~]=size(backbone);
mb=floor(m/spr);
if m<20; return; end
x=1:m;
xb=linspace(1,m,mb);
xb=xb';
x=x';
for k=1:3
    y=backbone(:,k);
    xd=[1,m,1,m];
    d=[y(1),y(m),y(2)-y(1),y(m)-y(m-1)];
    js=[0,0,1,1];
    pp=spfit(x,y,xb,3,xd,d,js);
    smoothed=ppval(pp,xb);
    if ~sum(isnan(smoothed))
        backbone(:,k)=interp1(xb,smoothed,x,'spline');
    end
end

function M = cardinal_spline(P,c,N)

s = (1-c)./2;

% Construct cardinal spline basis matrix
MC = [-s     2-s   s-2    s;
      2*s   s-3   3-2*s  -s;
      -s     0     s      0;
      0      1     0      0];

u = linspace(0,1,N+1).';
U = [u.^3 u.^2 u ones(size(u))];
M = U*MC*P;
