function [OE_osc] = me2osc(OE_m,J2,Re)
%me2osc Calculates the osculating orbital elements of a satellite
%   Based on algorithm in appendix F of "Analytical Mechanics of Space
%   Systems" Schaub & Junkins- Second Edition
%   Accepts a 6xN Matrix of osculating elements in the following order:
%
% a  - semimajor axis (km)
% e  - eccentricity
% i  - inclination relative to equatorial plane (deg)
% O  - right ascension of ascending node (deg)
% w  - argument of periapsis (deg)
% Me - Mean anomaly (deg)
% Outputs a 6xN Matrix on mean elements in the same order

if nargin == 1
    primary = earth();
    J2 = primary.J2;
    Re = primary.Re;
end

a0 = OE_m(1,:);
e0 = OE_m(2,:);
i0 = pi/180*OE_m(3,:);
O0 = wrapTo2Pi(pi/180*OE_m(4,:));
w0 = wrapTo2Pi(pi/180*OE_m(5,:));
M0 = wrapTo360(OE_m(6,:));
TA = me2ta(M0,e0);
M0 = pi/180*M0;
TA = pi/180*TA;

g2 = J2/2*(Re./a0).^2;

% Check inclination singularity
epsilon = 0.005;
den = 1-5*cos(i0).^2;
den = den.*(abs(den) >= epsilon) + epsilon*sign(den).*(abs(den) < epsilon);

% intermediate parameters
eta = sqrt(1-e0.^2);
g22 = g2./eta.^4;
a_r = (1+e0.*cos(TA))./eta.^2;
de1 = g22/8.*e0.*eta.^2.*(1-11*cos(i0).^2 - 40*cos(i0).^4./den).*cos(2*w0);

de = de1 + eta.^2/2.*(g2.*((3*cos(i0).^2 - 1)./eta.^6.*(e0.*eta + e0./(1+eta) + 3*cos(TA) ...
    + 3*e0.*cos(TA).^2 + e0.^2.*cos(TA).^3) + 3*(1 - cos(i0).^2)./eta.^6.*(e0 + 3*cos(TA) ...
    + 3*e0.*cos(TA).^2 + e0.^2.*cos(TA).^3).*cos(2*w0+2*TA)) ...
    - g22.*(1-cos(i0).^2).*(3*cos(2*w0 + TA) + cos(2*w0+3*TA)));

di = -e0.*de1./(eta.^2.*tan(i0)) + g22/2.*cos(i0).*sqrt(1-cos(i0).^2).*(3*cos(2*w0 + 2*TA) ...
    + 3*e0.*cos(2*w0+TA) + e0.*cos(2*w0 + 3*TA));

argsum = M0 + w0 + O0 ...
    + g22/8.*eta.^3.*(1 - 11*cos(i0).^2 - 40*cos(i0).^4./den) ...
    - g22/16.*(2 + e0.^2 - 11.*(2+3*e0.^2).*cos(i0).^2 ...
    - 40*(2 + 5*e0.^2).*cos(i0).^4./den - 400*e0.^2.*cos(i0).^6./den.^2)...
    + g22/4.*(-6*(1-5*cos(i0).^2).*(TA - M0 - e0.*sin(TA)) ...
    + (3-5*cos(i0).^2).*(3*sin(2*w0 + 2*TA) + 3*e0.*sin(2*w0+TA) ...
    + e0.*sin(2*w0 + 3*TA)))...
    - g22/8.*e0.^2.*cos(i0).*(11 + 80.*cos(i0).^2./den + 200*cos(i0).^4./den.^2)...
    - g22/2.*cos(i0).*(6*(TA - M0 + e0.*sin(TA)) ...
    - 3*sin(2*w0 + 2*TA) - 3*e0.*sin(2*w0+TA) - e0.*sin(2*w0 + 3*TA));

edM = g22/8.*e0.*eta.^3.*(1 - 11*cos(i0).^2 - 40*cos(i0).^4./den)...
    - g22/4.*eta.^3.*(2*(3*cos(i0).^2 - 1).*((a_r.*eta).^2 + a_r + 1).*sin(TA) ...
    + 3*(1-cos(i0).^2).*((-(a_r.*eta).^2 - a_r + 1).*sin(2*w0 + TA) ...
    + ((a_r.*eta).^2 + a_r + 1/3).*sin(2*w0 +3*TA)));

dO = -g22/8.*e0.^2.*cos(i0).*(11 + 80*cos(i0).^2./den + 200*cos(i0).^4./den.^2)...
    - g22/2.*cos(i0).*(6*(TA - M0 + e0.*sin(TA)) - 3.*sin(2*w0 + 2*TA)...
    - 3.*e0.*sin(2*w0 + TA) - e0.*sin(2*w0 + 3*TA));

% Compute Elements
d1 = (e0 + de).*sin(M0) + edM.*cos(M0);
d2 = (e0 + de).*cos(M0) - edM.*sin(M0);
d3 = (sin(i0/2) + cos(i0/2).*di/2).*sin(O0) + sin(i0/2).*dO.*cos(O0);
d4 = (sin(i0/2) + cos(i0/2).*di/2).*cos(O0) - sin(i0/2).*dO.*sin(O0);

a1 = a0 + a0.*g2.*((3.*cos(i0).^2-1).*(a_r.^3 - 1./eta.^3) ...
        + 3*(1-cos(i0).^2).*a_r.^3.*cos(2*w0 + 2*TA));

M1 = atan2(d1,d2);

e1 = sqrt(d1.^2 + d2.^2);

O1 = atan2(d3,d4);
% if any(sqrt(d3.^2 + d4.^2)>1)
%     bla = 4;
% end
i1 = 2*asin(sqrt(d3.^2 + d4.^2));

w1 = argsum - M1 - O1;

OE_osc = [a1; e1; wrapTo360(180/pi*i1);...
    wrapTo360(180/pi*O1); wrapTo360(180/pi*w1); wrapTo360(180/pi*M1)];
end