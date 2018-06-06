function x = get_RGT_orbit(a,e,inc,j,k,mu,J2,Rb,T_day)
%get_RGT_orbit derives one missing parameter for an RGT orbit
%   This function gets 2 orbital elements as well as the number of repeats:
%   j repeats in k days.
%   Calculates the remaining orbital element
%   The calculation takes into account only the effect of J2 perturbation
%   In order to choose which element to calculate, simply leave that field
%   as empty []
%   ~~~~~~~~~~~~~  INPUTS  ~~~~~~~~~~~~~~
%   a     - semimajor axis [km]
%   e     - eccentricity
%   inc   - orbital inclination [deg]
%   j     - number of repeats
%   k     - number of days to repeat
%   mu    - primary body gravitaional constant [km^3/s^2]
%   J2    - first zonal harmonic of primary
%   Rb    - primary body equatorial radius [km]
%   T_day - length of primary day [sec]

if nargin < 6 % default is Earth
    mu = 398600.440;
    J2 = 0.0010826265;
    Rb = 6378.137;
    T_day = 2*pi/7.29211585530e-5;
end
wb = 2*pi/T_day; % primary angular velocity
if isempty(a)
    a_init = (((k*T_day/j)/2/pi)^2*mu)^(1/3);
    syms a_var real;
    RAANdot = -3/2*J2*sqrt(mu)*Rb^2*a_var^(-7/2)*(1-e^2)^-2*cosd(inc);
    ArgPdot = 3/4*J2*sqrt(mu)*Rb^2*a_var^(-7/2)*(1-e^2)^-2*(5*(cosd(inc))^2 -1);
    MeAndot = 3/4*J2*sqrt(mu)*Rb^2*a_var^(-7/2)*(1-e^2)^(-3/2)*(3*(cosd(inc))^2 -1);
    eq = mu/a_var^3 == (j/k*(wb - RAANdot) - (MeAndot + ArgPdot))^2;
    x = double(vpasolve(eq,a_var,[Rb, 2*a_init]));
elseif isempty(e)
    syms e_var real;
    RAANdot = -3/2*J2*sqrt(mu)*Rb^2*a^(-7/2)*(1-e_var^2)^-2*cosd(inc);
    ArgPdot = 3/4*J2*sqrt(mu)*Rb^2*a^(-7/2)*(1-e_var^2)^-2*(5*(cosd(inc))^2 -1);
    MeAndot = 3/4*J2*sqrt(mu)*Rb^2*a^(-7/2)*(1-e_var^2)^(-3/2)*(3*(cosd(inc))^2 -1);
    eq = sqrt(mu)/a^(3/2) == j/k*(wb - RAANdot) - (MeAndot + ArgPdot);
    x = double(vpasolve(eq,e_var,[0,1]));
elseif isempty(inc)
    syms inc_var real;
    RAANdot = -3/2*J2*sqrt(mu)*Rb^2*a^(-7/2)*(1-e^2)^-2*cos(inc_var);
    ArgPdot = 3/4*J2*sqrt(mu)*Rb^2*a^(-7/2)*(1-e^2)^-2*(5*(cos(inc_var))^2 -1);
    MeAndot = 3/4*J2*sqrt(mu)*Rb^2*a^(-7/2)*(1-e^2)^(-3/2)*(3*(cos(inc_var))^2 -1);
    eq = sqrt(mu)/a^(3/2) == j/k*(wb - RAANdot) - (MeAndot + ArgPdot);
    x = 180/pi*double(vpasolve(eq,inc_var,[0,pi]));
else
    error('Orbit is already defined')
end
if isempty(x)
    error('No viable solution found')
end
end

