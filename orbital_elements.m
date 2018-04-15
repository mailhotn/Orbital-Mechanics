function [ e, h, i, O, w, th ] = orbital_elements( R, V ,mu )
%orbital_elements calculate orbital elements from state in geocentric
%equatorial reference frame
%
%  ~~~~~~~~~~   Variables    ~~~~~~~~~~~
% mu - gravitational parameter (km^3/s^2)
% e  - eccentricity
% h  - specific angular momentum (km^2/s)
% i  - inclination relative to equatorial plane (deg)
% O  - right ascension of ascending node (deg)
% w  - argument of periapsis (deg)
% th - true anomaly (deg)
if nargin < 3
    mu = 398600; % default is earth orbit of satellite with negligible mass
end
r = norm(R);
vr = dot(R,V)/r;
H = cross(R,V);
h = norm(H);
E = (cross(V,H)-mu*R/r)/mu;
e = norm(E);
i = acosd(H(3)/h);
N = cross([0 0 1],H);
n = norm(N);
if N(2) >= 0
    O = acosd(N(1)/n);
else
    O = 360 - acosd(N(1)/n);
end
if E(3) >= 0
    w = acosd(dot(N,E)/(n*e));
else
    w = 360 - acosd(dot(N,E)/(n*e));
end
if vr >= 0
    th = acosd(dot(E,R)/(e*r));
else
    th = 360 - acosd(dot(E,R)/(e*r));
end
end