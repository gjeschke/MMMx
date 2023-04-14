function lambda = debye_length(I,epsilon,T)
% lambda = debye_length(I)
%
% computes the Debye length at given ionic strength I (in mol/L)
%
% solvent permittivity espsilon can be provided or defaults to 78.2 (water)
% temperature can be given or defaults to 298 K
%
% output is in Angstrom
%
% G. Jeschke, 14.04.2023

if ~exist('epsilon','var') || isempty(epsilon)
    epsilon = 78.2;
end

if ~exist('T','var') || isempty(T)
    T = 298;
end

eps0 = 8.8542e-12;
kB = 1.380649e-23;
elcharge = 1.602176634e-19;
NA = 6.02214076e23;
I = 1000*I; % mol/L -> mol/m^3
lambda = 1e10*sqrt(eps0*epsilon*kB*T/(2*elcharge^2*NA*I));
