function [ a, e, i, O, w, th ] = state2elements( R, V ,mu )
%state2elements calculate orbital elements from state in geocentric
%equatorial reference frame. Accepts R and V matrices of N different points
%and outputs the elements as 1xN vectors
%
%  ~~~~~~~~~~   Variables    ~~~~~~~~~~~
% mu - gravitational parameter (km^3/s^2)
% e  - eccentricity
% a  - semimajor axis (km)
% i  - inclination relative to equatorial plane (deg)
% O  - right ascension of ascending node (deg)
% w  - argument of periapsis (deg)
% th - true anomaly (deg)
if nargin < 3
    mu = 398600; % default is earth orbit of satellite with negligible mass
end
[rows, ~] = size(R);
if rows ~= 3
    R = R.'; % Make sure R is a wide matrix
end
[rows, ~] = size(V);
if rows ~= 3
    V = V.'; % Make sure V is a wide matrix
end
r = sqrt(dot(R,R,1));
vr = dot(R,V,1)./r;
H = cross(R,V);
h = sqrt(dot(H,H,1));
E = (cross(V,H)-mu*R./r)/mu;
e = sqrt(dot(E,E,1));
i = acosd(H(3,:)./h);
N = cross(repmat([0 0 1].',1,size(H,2)),H);
n = sqrt(dot(N,N,1));

O = (N(2,:)>=0).*acosd(N(1,:)./n) + ...
    (N(2,:)<0).*(360 - acosd(N(1,:)./n));

w = (E(3,:)>=0).*acosd(dot(N,E,1)./(n.*e)) + ...
    (E(3,:)<0).*(360 - acosd(dot(N,E,1)./(n.*e)));

th = (vr>=0).*acosd(dot(E,R,1)./(e.*r)) + ...
    (vr<0).*(360 - acosd(dot(E,R,1)./(e.*r)));

a = h.^2./(mu*(1-e.^2));
end