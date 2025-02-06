function [xRsw] = eci2rsw(xEci,oe,primary)
%eci2rsw transforms coordinates from eci to rsw/lvlh, according to the orbital
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
if nargin < 3
    primary = earth();
end
if nargin < 2
    oe = eci2oe(xEci,[],primary,'me');
end

% a  = oe(1,:);
e  = oe(2,:);
i  = oe(3,:);
O  = oe(4,:);
w  = oe(5,:);
f = me2ta(oe(6,:),e);

xRsw = nan(size(xEci));
x = xEci(1,:);
y = xEci(2,:);
z = xEci(3,:);
dx = xEci(4,:);
dy = xEci(5,:);
dz = xEci(6,:);
% Perifocal
xPqw(1,:) = (cosd(O).*cosd(w)-sind(O).*sind(w).*cosd(i)).*x +...
    (sind(O).*cosd(w) + cosd(O).*cosd(i).*sind(w)).*y + ...
    sind(i).*sind(w).*z;
xPqw(2,:) = (-cosd(O).*sind(w) - sind(O).*cosd(i).*cosd(w)).*x +...
    (cosd(O).*cosd(i).*cosd(w) - sind(O).*sind(w)).*y +...
    cosd(w).*sind(i).*z;
xPqw(3,:) = sind(O).*sind(i).*x -cosd(O).*sind(i).*y + cosd(i).*z;
% velocities
xPqw(4,:) = (cosd(O).*cosd(w)-sind(O).*sind(w).*cosd(i)).*dx +...
    (sind(O).*cosd(w) + cosd(O).*cosd(i).*sind(w)).*dy + ...
    sind(i).*sind(w).*dz;
xPqw(5,:) = (-cosd(O).*sind(w) - sind(O).*cosd(i).*cosd(w)).*dx +...
    (cosd(O).*cosd(i).*cosd(w) - sind(O).*sind(w)).*dy +...
    cosd(w).*sind(i).*dz;
xPqw(6,:) = sind(O).*sind(i).*dx -cosd(O).*sind(i).*dy + cosd(i).*dz; 

% LVLH
xRsw(1,:) = cos(f).*xPqw(1,:) + sin(f).*xPqw(2,:);
xRsw(2,:) = -sin(f).*xPqw(1,:) + cos(f).*xPqw(2,:);
xRsw(3,:) = xPqw(3,:);

xRsw(4,:) = cos(f).*xPqw(4,:) + sin(f).*xPqw(5,:);
xRsw(5,:) = -sin(f).*xPqw(4,:) + cos(f).*xPqw(5,:);
xRsw(6,:) = xPqw(6,:);