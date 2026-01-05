classdef Moon
%Moon defines a class containing Moon physical characteristics
% all values are taken from Vallado
%
%  ~~~~~~~~~~   Primary Attributes    ~~~~~~~~~~~
% mu - gravitational parameter (km^3/s^2)
% Re - equatorial radius (km)
% J2 - zonal harmonic
% we - rotation rate (deg/s)
% Tday - length of sidereal day (s)
% Tyear - length of year (s) - orbital period
%
%  ~~~~~~~~~~   Third-Body Attributes    ~~~~~~~~~~~
% sma - semimajor axis (km)
% ecc - eccentricity
% inc - inclination (deg)
% nMo - mean motion (1/s)
properties (Constant)
name = 'Moon';

% Primary Attributes
mu = 4902.799;
Re = 1738.0;
R = 1738.0;
J2 = 2.027e-4;
we = 360/27.32166/86400;
Tday = 27.32166*86400;
Tyear = 27.321582*86400;

% Third-body Atributes (earth primary)
sma = 384400;
ecc = 0.05490;
inc = 5.145396;
nMo = 2*pi/27.321582/86400;
mass = 7.3483e22;



end
end