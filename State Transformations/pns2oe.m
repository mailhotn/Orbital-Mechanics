function [oe] = pns2oe (pns, primary)
%pns2oe transforms from NONSINGULAR polar nodal to classical elements
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
% th - modified argument of latitude (rad)
% n  - modified right ascension of the ascending node (rad)
% R  - radial velocity (km/s)
% TH - specific angular momentum (km^2/s^2)
% N  - modified angular momentum in Z axis
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

if nargin < 2
    primary = earth;
end
mu = primary.mu;
% check input size
if size(pns,1) ~=6 
    if size(pns,2) == 6
        pns = pns.';
    else
        error('Bad input - no dimension is 6')
    end
end

r = pns(1,:);
th = pns(2,:);
n = pns(3,:);
R = pns(4,:);
TH = pns(5,:);
N = pns(6,:);

k = TH.^2./(mu*r)-1;
s = TH.*R/mu;
f = 180/pi*unwrap(atan2(s,k));
% remove singularities n/N = nan only if N=0 i.e. orbit is equatorial
nbyN = n./N;
nbyN(isnan(nbyN)) = 0;

oe(1,:) = mu*r.^2./(2*mu*r -R.^2.*r.^2 -TH.^2);
oe(2,:) = sqrt(TH.^4 +r.*(R.^2.*r -2*mu).*TH.^2 +mu^2*r.^2)./(mu*r);
oe(3,:) = acosd(1-N.^2./(2*TH));
oe(4,:) = wrapTo360(-180/pi*nbyN);
oe(5,:) = wrapTo360(180/pi*th - f + 180/pi*nbyN);
oe(6,:) = wrapTo360(f);


