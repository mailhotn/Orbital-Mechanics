function earthstruct = earth()
%earth Returns a struct containing Earth physical characteristics
% all values are taken from JPL HORIZONS database
%
%  ~~~~~~~~~~   Attributes    ~~~~~~~~~~~
% mu - gravitational parameter (km^3/s^2)
% R  - equatorial radius (km)
% J2 - zonal harmonic

earthstruct.mu = 398600.440;
earthstruct.R  = 6378.137;
earthstruct.J2 = 0.0010826265;

end

