function OE = eci2oe( X1, X2 ,primary )
%state2elements calculate orbital elements from state in geocentric
%equatorial reference frame. Accepts R and V matrices of N different points
%and outputs the elements as a 6xN matrix
%
% Also accepts full 6xN state matrix
% outputs:    OE = [a; e; i; O; w; th];
%  ~~~~~~~~~~   Variables    ~~~~~~~~~~~
% mu - gravitational parameter (km^3/s^2)
% e  - eccentricity
% a  - semimajor axis (km)
% i  - inclination relative to equatorial plane (deg)
% O  - right ascension of ascending node (deg)
% w  - argument of periapsis (deg)
% th - true anomaly (deg)
if nargin < 2
    X2 = [];
end
if nargin < 3
    primary = earth();
     % default is earth orbit of satellite with negligible mass
end
mu = primary.mu;
if ~isempty(X2) % X1 is R, X2 is V
    R = X1;
    V = X2;
else % X2 not given, X1 is full state
    if size(X1,1) ~= 6 && size(X1,2) == 6
        X1 = X1.';
    end
    R = X1(1:3,:);
    V = X1(4:6,:);
end
if size(R,1) ~= 3
    R = R.'; % Make sure R is a wide matrix
end
if size(V,1) ~= 3
    V = V.'; % Make sure V is a wide matrix
end
r = vecnorm(R,2);
vr = dot(R,V,1)./r;
H = cross(R,V);
h = vecnorm(H,2);
E = (cross(V,H)-mu*R./r)/mu;
e = vecnorm(E,2);
i = acosd(H(3,:)./h);
N = cross(repmat([0 0 1].',1,size(H,2)),H);
n = vecnorm(N,2);

O = (N(2,:)>=0).*acosd(N(1,:)./n) + ...
    (N(2,:)<0).*(360 - acosd(N(1,:)./n));

w = (E(3,:)>=0).*acosd(dot(N./n,E./e,1)) + ...
    (E(3,:)<0).*(360 - acosd(dot(N./n,E./e,1)));

th = (vr>=0).*acosd(dot(E./e,R./r,1)) + ...
    (vr<0).*(360 - acosd(dot(E./e,R./r,1)));

a = h.^2./(mu*(1-e.^2));
OE = real([a; e; i; O; w; th]);
% Just get rid of imaginary stuff. For some reason matlab gets and error of
% O(eps) when calculating the product of two normalized vectors, there is
% no reason to get imaginary numbers from the acos
end