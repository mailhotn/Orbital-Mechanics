function oe3b = eci2oe3b(rEci3b,vEci3b,third)
% converts a position vector of third body to orbital elements

sma = third.sma;
ecc = third.ecc;

H = cross(rEci3b,vEci3b);
h = vecnorm(H,2);
inc = acosd(H(3,:)./h);
N = cross(repmat([0 0 1].',1,size(H,2)),H);
n = vecnorm(N,2);

ran = (N(2,:)>=0).*acosd(N(1,:)./n) + ...
    (N(2,:)<0).*(360 - acosd(N(1,:)./n));

aop = (E(3,:)>=0).*acosd(dot(N./n,E./e,1)) + ...
    (E(3,:)<0).*(360 - acosd(dot(N./n,E./e,1)));

tan = (vr>=0).*acosd(dot(E./e,R./r,1)) + ...
    (vr<0).*(360 - acosd(dot(E./e,R./r,1)));

oe3b = [sma;ecc;inc;ran;aop;tan];