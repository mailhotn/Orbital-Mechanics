function earthstruct = earth()
%earth Returns a struct containing Earth physical characteristics
% all values are taken from JPL HORIZONS database
% ~~~~~~~~~~~ Might be better as a Class ~~~~~~~~~~~
%
%  ~~~~~~~~~~   Attributes    ~~~~~~~~~~~
% mu - gravitational parameter (km^3/s^2)
% R  - equatorial radius (km)
% J2 - zonal harmonic
% we - earth rotation rate (deg/s)
% Tday - length of sidereal day (s)

earthstruct.mu = 398600.440;
earthstruct.Re = 6378.137;
earthstruct.J2 = 0.0010826265;
earthstruct.we = 180/pi*7.29211585530e-5;
earthstruct.Tday = 2*pi/7.29211585530e-5;

end

