function xEci = rsw2eci(xRsw,oe)
%rsw2eci transforms coordinates from rsw/lvlh to eci, according to the orbital
% elements defining the rotation.
% input is a 6xN matrix of orbital elements:
%
%        oe = [a; e; i; O; w; M];
%
%  ~~~~~~~~~~   Variables   ~~~~~~~~~~
% xEci - 6xN matrix of ECI vectors
% e  - eccentricity
% a  - semi-major axis (km)
% i  - inclination relative to equatorial plane (deg)
% O  - right ascension of ascending node (deg)
% w  - argument of periapsis (deg)
% M  -  mean anomaly (deg)
%
%  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%


% a  = oe(1,:);
ecc  = oe(2,:);
inc  = oe(3,:);
ran  = oe(4,:);
aop  = oe(5,:);
f = me2ta(oe(6,:),ecc);

xEci = nan(size(xRsw));
r = xRsw(1,:);
s = xRsw(2,:);
w = xRsw(3,:);
dr = xRsw(4,:);
ds = xRsw(5,:);
dw = xRsw(6,:);

% LVLH to perifocal
p = cosd(f).*r -sind(f).*s;
q = sind(f).*r + cosd(f).*s;

dp = cosd(f).*dr -sind(f).*ds;
dq = sind(f).*dr + cosd(f).*ds;

% Perifocal to eci
xEci(1,:) = (cosd(ran).*cosd(aop)-sind(ran).*sind(aop).*cosd(inc)).*p +...
    (-cosd(ran).*sind(aop) - sind(ran).*cosd(aop).*cosd(inc)).*q + ...
    sind(inc).*sind(ran).*w;
xEci(2,:) = (sind(ran).*cosd(aop) + cosd(ran).*cosd(inc).*sind(aop)).*p +...
    (cosd(ran).*cosd(inc).*cosd(aop) - sind(ran).*sind(aop)).*q +...
    -cosd(ran).*sind(inc).*w;
xEci(3,:) = sind(aop).*sind(inc).*p +cosd(aop).*sind(inc).*q + cosd(inc).*w;
% velocities
xEci(4,:) = (cosd(ran).*cosd(aop)-sind(ran).*sind(aop).*cosd(inc)).*dp +...
    (-cosd(ran).*sind(aop) - sind(ran).*cosd(aop).*cosd(inc)).*dq + ...
    sind(inc).*sind(ran).*dw;
xEci(5,:) = (sind(ran).*cosd(aop) + cosd(ran).*cosd(inc).*sind(aop)).*dp +...
    (cosd(ran).*cosd(inc).*cosd(aop) - sind(ran).*sind(aop)).*dq +...
    -cosd(ran).*sind(inc).*dw;
xEci(6,:) = sind(aop).*sind(inc).*dp +cosd(aop).*sind(inc).*dq + cosd(inc).*dw;

