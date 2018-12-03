function earthstruct = normalized()
%earth Returns a struct containing physical characteristics of a normalized
% planet
% ~~~~~~~~~~~ Might be better as a Class ~~~~~~~~~~~
%
%  ~~~~~~~~~~   Attributes    ~~~~~~~~~~~
% mu - gravitational parameter (km^3/s^2)
% R  - equatorial radius (km)
% J2 - zonal harmonic
% we - earth rotation rate (deg/s)
% Tday - length of sidereal day (s)

earthstruct.mu = 1;
earthstruct.Re = 1;
end