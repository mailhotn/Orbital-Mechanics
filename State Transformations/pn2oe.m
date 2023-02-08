function [oe] = pn2oe (pn, primary)
%pn2oe transforms from polar nodal to classical elements
% angle inputs in radians, output in degrees
% input is a 6xN matrix of orbital elements:
%
%        OE = [a; e; i; O; w; f];
%        pn = [r; th; n; R; TH; N];
%  ~~~~~~~~~~   Variables   ~~~~~~~~~~
% mu - gravitational parameter (km^3/s^2)
% e  - eccentricity
% a  - semi-major axis (km)
% i  - inclination relative to equatorial plane (deg)
% O  - right ascension of ascending node (deg)
% w  - argument of periapsis (deg)
% f  - true anomaly (deg)
%
% r  - radial position (km)
% th - argument of latitude (rad)
% n  - right ascension of the ascending node (rad)
% R  - radial velocity (km/s)
% TH - specific angular momentum (km^2/s^2)
% N  - angular momentum in Z axis (km^2/s^2)
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if nargin < 2
    primary = earth;
end
mu = primary.mu;
% check input size
if size(pn,1) ~=6 
    if size(pn,2) == 6
        pn = pn.';
    else
        error('Bad input - no dimension is 6')
    end
end

r = pn(1,:);
th = pn(2,:);
n = pn(3,:);
R = pn(4,:);
TH = pn(5,:);
N = pn(6,:);

k = TH.^2./(mu*r)-1;
s = TH.*R/mu;
f = 180/pi*unwrap(atan2(s,k));

oe(1,:) = mu*r.^2./(2*mu*r -R.^2.*r.^2 -TH.^2);
oe(2,:) = sqrt(TH.^4 +r.*(R.^2.*r -2*mu).*TH.^2 +mu^2*r.^2)./(mu*r);
oe(3,:) = acosd(N./TH);
oe(4,:) = 180/pi*n;
oe(5,:) = 180/pi*th - f;
oe(6,:) = f;


