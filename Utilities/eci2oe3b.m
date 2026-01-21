function [oe3b,E] = eci2oe3b(rEci3b,vEci3b,primary,third)
% converts a position vector of third body to orbital elements

% Make sure the matrices are wide
if size(rEci3b,1) ~=3
    rEci3b = rEci3b.';
end
if size(vEci3b,1) ~=3
    vEci3b = vEci3b.';
end
nT = size(rEci3b,2);

mu = primary.mu+third.mu;
sma = third.sma;
ecc = third.ecc;

r = vecnorm(rEci3b,2);
vr = dot(rEci3b ...
    ,vEci3b,1)./r;
H = cross(rEci3b,vEci3b);
h = vecnorm(H,2);
E = (cross(vEci3b,H)-mu*rEci3b./r)/mu;
ecc = vecnorm(E,2);
inc = acos(H(3,:)./h);
N = cross(repmat([0 0 1].',1,size(H,2)),H);
n = vecnorm(N,2);

ran = (N(2,:)>=0).*acos(N(1,:)./n) + ...
    (N(2,:)<0).*(2*pi - acos(N(1,:)./n));

aop = (E(3,:)>=0).*acos(dot(N./n,E./ecc,1)) + ...
    (E(3,:)<0).*(2*pi - acos(dot(N./n,E./ecc,1)));

tan = (vr>=0).*acos(dot(E./ecc,rEci3b./r,1)) + ...
    (vr<0).*(2*pi - acos(dot(E./ecc,rEci3b./r,1)));

sma = h.^2./(mu*(1-ecc.^2));
% sma = third.sma*ones(1,nT);
% ecc = third.ecc*ones(1,nT);


oe3b = [sma;ecc;inc;ran;aop;tan];