function sunStruct = sun()
%sun Returns a struct containing Earth physical characteristics
% all values are taken from Vallado
% ~~~~~~~~~~~ Might be better as a Class ~~~~~~~~~~~
%
%  ~~~~~~~~~~   Attributes    ~~~~~~~~~~~
% mu - gravitational parameter (km^3/s^2)
% Re - equatorial radius (km)
% J2 - zonal harmonic
% we - earth rotation rate (deg/s)
% Tday - length of sidereal day (s)
% Tyear - length of year (s)
% ExpAtmoModel - [h_b h_c rho_b(kg/m^3)
%
%  ~~~~~~~~~~   Third-Body Attributes    ~~~~~~~~~~~
% sma - semimajor axis (km)
% ecc - eccentricity
% inc - inclination (deg)
% nMo - mean motion (1/s)

sunStruct.mu = 398600.4415;
sunStruct.Re = 6378.1363;
sunStruct.J2 = 1.08262617385222e-3;
sunStruct.we = 180/pi*7.29211585530e-5;
sunStruct.Tday = 2*pi/7.29211585530e-5;
sunStruct.Tyear = 365.2421897*24*60*60;


% Third-body Atributes (moon primary)
sunStruct.sma = 149598023.0;
sunStruct.ecc = 0.05490;
sunStruct.inc = 5.145396+6.68;
sunStruct.nMo = 2*pi/sunStruct.Tyear;
sunStruct.mass = 1.9891e30;

end