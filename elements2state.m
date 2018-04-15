function [R, V] = elements2state (a, e, i, O, w, th, mu)
%elements2state calculates Position & velocity vectors in ECI Frame
% from the orbital elements
%
%  ~~~~~~~~~~   Variables   ~~~~~~~~~~
% mu - gravitational parameter (km^3/s^2)
% e  - eccentricity
% a  - semi-major axis (km)
% i  - inclination relative to equatorial plane (deg)
% O  - right ascension of ascending node (deg)
% w  - argument of periapsis (deg)
% th - true anomaly (deg)
%
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% Note: th can be a vector in which case R and V will be 3 x N wide
% matrices. This feature is useful for drawing orbit segments.
% 
if nargin < 7
    mu = 398600; % default is earth orbit of satellite with negligible mass
end
[rows, cols] = size(th);
if rows > cols
    th = th.'; % Make sure th is a row vector
end
r = a*(1-e^2)./(1+e*cosd(th)).*[cosd(th); sind(th); zeros(1,length(th))];
v_p = sqrt(mu/a/(1-e^2)).*(1+e*cosd(th));
v_r = sqrt(mu/a/(1-e^2)).*(e*sind(th));
u_p = [-sind(th); cosd(th); zeros(1,length(th))];
u_r = [cosd(th); sind(th); zeros(1,length(th))];
v = v_p.*u_p +v_r.*u_r;
Rot = [ cosd(O)*cosd(w) - sind(O)*cosd(i)*sind(w), - cosd(O)*sind(w) - sind(O)*cosd(i)*cosd(w),  sind(O)*sind(i);
        sind(O)*cosd(w) + cosd(O)*cosd(i)*sind(w),   cosd(O)*cosd(i)*cosd(w) - sind(O)*sind(w), -cosd(O)*sind(i);
                                  sind(i)*sind(w),                             cosd(w)*sind(i),         cosd(i)];                      
R = Rot*r;  
V = Rot*v;
end