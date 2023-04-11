function coefficients = point_spread(element)
% coefficients = point_spread(element)
%
% returns set of four pairs of Gaussian point-spread coefficients for
% electron density maps for an element specified by its atom number
%
% the coefficients are as defined in the SI of:
% R. Briones, C. Blau, C. Kutzner, B. L. de Groot, C. Aponte-Santamarı,
% Biophys. J. 116, 4–11 (2019). https://doi.org/10.1016/j.bpj.2018.11.3126
%
% data were taken from file electronscattering.dat in the GROMACS extension
% available at https://mptg-cbp.github.io/gromaps.html
%
% G. Jeschke, 11.04.2023

H = [3.360350e-02,8.435271e+01;9.360282e-02,2.468311e+02;1.027241e-01,2.482872e+02;3.532598e-01,1.697996e+03];    
C = [9.253835e-02,6.890005e+01;6.260279e-01,1.834891e+02;1.048004e+00,6.026779e+02;3.033389e-01,2.226268e+03];
N = [8.946127e-02,9.246527e+01;5.702793e-01,2.362661e+02;1.181245e+00,7.689153e+02;8.856041e-01,1.112068e+04];
O = [2.126579e-02,8.683121e+01;2.289285e-01,1.794531e+02;1.006541e+00,5.041622e+02;2.285777e+00,4.493332e+03];
P = [1.492938e-01,5.156540e+01;9.585026e-01,1.223678e+02;1.564325e+00,3.209862e+02;1.646213e+00,1.572971e+03];
S = [2.272068e-01,7.166818e+01;1.373602e+00,1.850285e+02;1.905039e+00,6.119167e+02;6.346896e-01,2.793350e+03];

switch element
    case 1
        coefficients = H;
    case 6
        coefficients = C;
    case 7
        coefficients = N;
    case 8
        coefficients = O;
    case 15
        coefficients = P;
    case 16
        coefficients = S;
    otherwise
        coefficients = C;
end
% conversion of B coefficients from nm^-2 to Å^-2
coefficients(:,2) = coefficients(:,2)/100;