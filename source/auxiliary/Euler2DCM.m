function dircos=Euler2DCM(Euler)
% Converts Euler angles in radians to a direction cosine matrix dircos to 
% Euler angles refer to rotations about the z, y', and z'' axis, as usual
% in magnetic resonance contexts
% function is just a wrapper for SpinCalc by John Fuller
%
% G. Jeschke, 2009

dircos=SpinCalc('EA323toDCM',180*Euler/pi,0.1,0);
