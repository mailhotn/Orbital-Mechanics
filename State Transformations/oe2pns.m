function [pns] = oe2pns (oe, primary)
%oe2pns transforms from classical to NONSINGULAR polar nodal elements
% angle inputs in degrees, output in radians
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
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
if nargin < 2
    primary = earth;
end
mu = primary.mu;
% check input size
if size(oe,1) ~=6 
    if size(oe,2) == 6
        oe = oe.';
    else
        error('Bad input - no dimension is 6')
    end
end
sma = oe(1,:);
ecc = oe(2,:);
inc = pi/180*oe(3,:);
ran = pi/180*oe(4,:);
aop = pi/180*oe(5,:);
f = pi/180*oe(6,:);

p = sma.*(1-ecc.^2);
N = sqrt(mu*p).*cos(inc);
h = sqrt(mu*p);
pns = nan(size(oe));

pns(1,:) = p./(1+ecc.*cos(f));
pns(2,:) = f + aop + ran;
pns(3,:) = -sqrt(2*h-2*N).*ran;
pns(4,:) = mu*ecc.*sin(f)./sqrt(mu*p);
pns(5,:) = h;
pns(6,:) = sqrt(2*h-2*N);


