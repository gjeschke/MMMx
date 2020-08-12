function transmat=Euler_affine_fast(angvec,transmat)
% function transmat=Euler_affine_fast(angvec,mat)
%
% Creates a 4x4 matrix for an affine coordinate transformation in 3D space
% see: M. Bender, M. Brill, Computergrafik, Hanser, München, 2. Aufl.,
%      2006, section 2.1
% fast version for Monte Carlo trials in helix bundle assembler, created
% from affine.m
% angvec           Euler angles [alpha,beta,gamma]
%                  for rotations about z, y', and
%                  z'' axes, in radians
%                  WARNING: the objects are
%                  rotated, not the frame, use
%                  -alpha, -beta, -gamma for the
%                  frame rotation more usual in
%                  magnetic resonance
% transmat         must be an eye(4) matrix, providing it obviates
%                  expensive initialization
%
% G. Jeschke, 2012

ca=cos(angvec(1)); sa=sin(-angvec(1));
cb=cos(angvec(2)); sb=sin(-angvec(2));
cg=cos(angvec(3)); sg=sin(-angvec(3));
cgcb=cg*cb;
msgcb=-sg*cb;
transmat(1,1)=cgcb*ca-sg*sa;
transmat(1,2)=cgcb*sa+sg*ca;
transmat(1,3)=-cg*sb;
transmat(2,1)=msgcb*ca-cg*sa;
transmat(2,2)=msgcb*sa+cg*ca;
transmat(2,3)=sg*sb;
transmat(3,1)=sb*ca;
transmat(3,2)=sb*sa;
transmat(3,3)=cb;
       