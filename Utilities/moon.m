function moonStruct = moon()
%moon Returns a struct containing Moon physical characteristics
% all values are taken from Vallado
% ~~~~~~~~~~~ Might be better as a Class ~~~~~~~~~~~
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

moonStruct.name = 'Moon';

% Primary Attributes
moonStruct.mu = 4902.799;
moonStruct.Re = 1738.0;
moonStruct.J2 = 2.027e-4;
moonStruct.we = 360/27.32166/86400;
moonStruct.Tday = 27.32166*86400;
moonStruct.Tyear = 27.321582*86400;

% Third-body Atributes (earth primary)
moonStruct.sma = 384400;
moonStruct.ecc = 0.05490;
moonStruct.inc = 5.145396;
moonStruct.nMo = 2*pi/moonStruct.Tyear;
moonStruct.mass = 7.3483e22;



end