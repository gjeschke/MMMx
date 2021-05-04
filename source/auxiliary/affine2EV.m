function EV=affine2EV(transmat)
% function EV=affine2EV(transmat)
%
% extracts the rotation matrix (direction cosine matrix) DCM from a 4x4
% affine transformation matrix and converts it to an Euler vector and
% corresponding Euler angle
% checks, whether the "rotation block" of the affine transformation matrix
% corresponds to only a rotation (set of three orthonormal vectors) with
% determinant 1
% otherwise, an empty matrix is returned
% this is just a wrapper for SpinCalc with sign correction
%
% G. Jeschke, 2010

tol=1e-3; % tolerance for the norm of the direction vectors
tol_proj=1e-3; % tolerance for cosine of angle between direction vectors
tol_det=1e-3; % tolerance for determinant

DCM=[];

DCM0=transmat(1:3,1:3);

if abs(det(DCM0)-1)>tol_det, return; end;

for k=1:3,
    if abs(norm(DCM0(k,:))-1) > tol, return; end;
    if abs(norm(DCM0(:,k))-1) > tol, return; end;
end;

for k=1:2,
    for kk=k+1:3,
        proj=abs(sum(DCM0(k,:).*DCM0(kk,:)));
        if proj>tol_proj; return; end;
        proj=abs(sum(DCM0(:,k).*DCM0(:,kk)));
        if proj>tol_proj; return; end;
    end;
end;

DCM=DCM0;

EV=SpinCalc('DCMtoEV',DCM,0.001,0);
EV(1:3)=-EV(1:3);