function marsstruct = mars()
%mars Returns a struct containing Mars physical characteristics
% all values are taken from Vallado
% ~~~~~~~~~~~ Might be better as a Class ~~~~~~~~~~~
%
%  ~~~~~~~~~~   Attributes    ~~~~~~~~~~~
% mu - gravitational parameter (km^3/s^2)
% R  - equatorial radius (km)
% J2 - zonal harmonic
% we - earth rotation rate (deg/s)
% Tday - length of sidereal day (s)
% Tyear - length of year (s)
% ExpAtmoModel - [h_b h_c rho_b(kg/m^3)

marsstruct.mu = 4.305e4;
marsstruct.Re = 3397.2;
marsstruct.J2 = 0.001964;
marsstruct.J3 = 0.000036;
marsstruct.Tday = 1.02395675*(24*60*60);
marsstruct.we = 360/(1.02395675*(24*60*60));
marsstruct.Tyear = 686.9150*24*60*60;
marsstruct.dOs = 2*pi/(686.9150*24*60*60);

end