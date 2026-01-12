classdef Sun
%Sun defines a class containing Sun physical characteristics
% all values are taken from Vallado
%
%  ~~~~~~~~~~   Attributes    ~~~~~~~~~~~
% mu - gravitational parameter (km^3/s^2)
% Re - equatorial radius (km)
% J2 - zonal harmonic
%
%  ~~~~~~~~~~   Third-Body Attributes    ~~~~~~~~~~~
% sma - semimajor axis (km)
% ecc - eccentricity
% inc - inclination (deg)
% nMo - mean motion (1/s)

properties (Constant)
name = "Sun";
mu = 1.32712428e11;
Re = 696000;
R = 696000;

% Third-body Attributes
sma = 149598023.0;
ecc = 0.016708617;
inc = 0;
nMo = 2*pi/(365.2421897*86400);
mass = 1.9891e30;

end
end