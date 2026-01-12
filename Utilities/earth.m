classdef Earth
%Earth defines a class containing Earth physical characteristics
% all values are taken from Vallado
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
properties (Constant)
name = 'Earth';

mu = 398600.4415;
Re = 6378.1363;
R = 6378.1363;
J2 = 1.08262617385222e-3;
J3 = -2.53241051856772e-6;
we = 180/pi*7.29211585530e-5;
Tday = 2*pi/7.29211585530e-5;
Tyear = 365.2421897*24*60*60;
dOs = 2*pi/(365.2421897*24*60*60);
ExpAtmoModel = ...
     [ 0   25 1.225e-00  7.249;
      25   30 3.899e-02  6.349;
      30   40 1.774e-02  6.682;
      40   50 3.972e-03  7.554;
      50   60 1.057e-03  8.382;
      60   70 3.206e-04  7.714;
      70   80 8.770e-05  6.549;
      80   90 1.905e-05  5.799;
      90  100 3.396e-06  5.382;
     100  110 5.297e-07  5.877;
     110  120 9.661e-08  7.263;
     120  130 2.438e-08  9.473;
     130  140 8.484e-09 12.636;
     140  150 3.845e-09 16.149;
     150  180 2.070e-09 22.523;
     180  200 5.464e-10 29.740;
     200  250 2.789e-10 37.105;
     250  300 7.248e-11 45.546;
     300  350 2.418e-11 53.628;
     350  400 9.518e-12 53.298;
     400  450 3.725e-12 58.515;
     450  500 1.585e-12 60.828;
     500  600 6.967e-13 63.822;
     600  700 1.454e-13 71.835;
     700  800 3.614e-14 88.667;
     800  900 1.170e-14 124.64;
     900 1000 5.245e-15 181.05;
    1000  inf 3.019e-15 268.00];

% Third-body Attributes (moon primary)
sma = 384400;
ecc = 0.05490;
inc = 5.145396+6.68;
nMo = 2*pi/(27.321582*86400);
mass = 5.9742e24;

end
end