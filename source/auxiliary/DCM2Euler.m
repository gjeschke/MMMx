function Euler=DCM2Euler(dircos)
% Converts a direction cosine matrix dircos to Euler angles in radians,
% Euler angles refer to rotations about the z, y', and z'' axis, as usual
% in magnetic resonance contexts
% function is just a wrapper for SpinCalc by John Fuller
%
% G. Jeschke, 2009

Euler=pi*SpinCalc('DCMtoEA323',dircos,0.1,0)/180;
for k=1:length(Euler),
    if Euler(k)>pi,
        Euler(k)=Euler(k)-2*pi;
    end;
end;