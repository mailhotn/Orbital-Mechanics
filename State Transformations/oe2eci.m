function [R, V] = oe2eci (OE, mu)
%elements2state calculates Position & velocity vectors in ECI Frame
% from the orbital elements
% input is a 6xN matrix of orbital elements:
%
%        OE = [a; e; i; O; w; th];
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
% 
if nargin < 2
    mu = 398600.440; % default is earth orbit of satellite with negligible mass
end
if size(OE,1) ~=6
    OE = OE.';
end
R = zeros(3,size(OE,2));
V = zeros(3,size(OE,2));
a  = OE(1,:);
e  = OE(2,:);
i  = OE(3,:);
O  = OE(4,:);
w  = OE(5,:);
th = OE(6,:);

r = a.*(1-e.^2)./(1+e.*cosd(th)).*[cosd(th); sind(th); zeros(1,length(th))];
v_p = sqrt(mu./a./(1-e.^2)).*(1+e.*cosd(th));
v_r = sqrt(mu./a./(1-e.^2)).*(e.*sind(th));
u_p = [-sind(th); cosd(th); zeros(1,length(th))];
u_r = [cosd(th); sind(th); zeros(1,length(th))];
v = v_p.*u_p +v_r.*u_r;

% rotate from pqw to ECI
R(1,:) = (cosd(O).*cosd(w) - sind(O).*cosd(i).*sind(w)).*r(1,:) + ...
         (-cosd(O).*sind(w) - sind(O).*cosd(i).*cosd(w)).*r(2,:) + ...
         (sind(O).*sind(i)).*r(3,:);
R(2,:) = (sind(O).*cosd(w) + cosd(O).*cosd(i).*sind(w)).*r(1,:) + ...
         (cosd(O).*cosd(i).*cosd(w) - sind(O).*sind(w)).*r(2,:) + ...
         (-cosd(O).*sind(i)).*r(3,:);
R(3,:) = (sind(i).*sind(w)).*r(1,:) + ...
         (cosd(w).*sind(i)).*r(2,:) + ...
         (cosd(i)).*r(3,:);
     
V(1,:) = (cosd(O).*cosd(w) - sind(O).*cosd(i).*sind(w)).*v(1,:) + ...
         (-cosd(O).*sind(w) - sind(O).*cosd(i).*cosd(w)).*v(2,:) + ...
         (sind(O).*sind(i)).*v(3,:);
V(2,:) = (sind(O).*cosd(w) + cosd(O).*cosd(i).*sind(w)).*v(1,:) + ...
         (cosd(O).*cosd(i).*cosd(w) - sind(O).*sind(w)).*v(2,:) + ...
         (-cosd(O).*sind(i)).*v(3,:);
V(3,:) = (sind(i).*sind(w)).*v(1,:) + ...
         (cosd(w).*sind(i)).*v(2,:) + ...
         (cosd(i)).*v(3,:);
% Previous rotation mechanism:
% Rot = [ cosd(O).*cosd(w) - sind(O).*cosd(i).*sind(w), - cosd(O).*sind(w) - sind(O).*cosd(i).*cosd(w),  sind(O).*sind(i);
%         sind(O).*cosd(w) + cosd(O).*cosd(i).*sind(w),   cosd(O).*cosd(i).*cosd(w) - sind(O).*sind(w), -cosd(O).*sind(i);
%                                   sind(i).*sind(w),                             cosd(w).*sind(i),         cosd(i)];                      
% R = Rot*r;  
% V = Rot*v;
end