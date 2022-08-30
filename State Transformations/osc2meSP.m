function [oeM] = osc2meSP(oeOsc,primary)
%osc2meSP Calculates the osculating orbital elements of a satellite
%   Based on Kozai 1959
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
end
    

J2 = primary.J2;
Re = primary.Re;

sma0 = oeOsc(1,:);
ecc0 = oeOsc(2,:);
inc0 = pi/180*oeOsc(3,:);
ran0 = wrapTo2Pi(pi/180*oeOsc(4,:));
aop0 = wrapTo2Pi(pi/180*oeOsc(5,:));
man0 = wrapTo360(oeOsc(6,:));
f = me2ta(man0,ecc0);
man0 = pi/180*man0;
f = pi/180*f;

A2 = -J2*Re^2*1.5; % positive A2: osc2me

% sma0 = sma0.*(1-A2./(sma0.*(1-ecc0.^2)).^2.*(1-1.5*sin(inc0).^2).*eta); % Idk

aR = (1+ecc0.*cos(f))./(1-ecc0.^2);
eta = sqrt(1-ecc0.^2);
p = sma0.*(1-ecc0.^2);
s2i = sin(inc0).^2;

% Short period variations
dSmaS = A2./sma0.*(2/3*(1-1.5*s2i).*(aR.^3 - 1./eta.^3) + ...
        aR.^3.*s2i.*cos(2*f+2*aop0));

dEccS = (1-ecc0.^2)./ecc0*A2./sma0.^2.*((1-1.5*s2i)/3.*(aR.^3-eta.^3) + ...
        0.5*aR.^3.*s2i.*cos(2*f+2*aop0)) - ...
        s2i/2./ecc0*A2./(sma0.*p).*(cos(2*f+2*aop0) +...
        ecc0.*cos(f+2*aop0) + ecc0/3.*cos(3*f+2*aop0));
    
dIncS = 0.25*A2./p.^2.*sin(2*inc0).*(cos(2*f+2*aop0) + ecc0.*cos(f+2*aop0) ...
    + ecc0/3.*cos(3*f+2*aop0));

dRanS = -A2./p.^2.*cos(inc0).*(f-man0+ecc0.*sin(f) -0.5*sin(2*f+2*aop0) ...
    -0.5*ecc0.*sin(f+2*aop0) -ecc0/6.*sin(3*f+2*aop0));
    
dAopS = A2./p.^2.*((2-2.5*s2i).*(f-man0+ecc0.*sin(f)) + ...
    (1-1.5*s2i).*((1./ecc0 - 0.25*ecc0).*sin(f) + 0.5*sin(2*f) + ecc0/12.*sin(3*f))...
    -(0.25./ecc0.*s2i + (0.5-15/16*s2i).*ecc0).*sin(f+2*aop0) ...
    +ecc0/16.*s2i.*sin(f-2*aop0) -0.5*(1-2.5*s2i).*sin(2*f+2*aop0)...
    +(7/12*s2i./ecc0-(1-19/8*s2i).*ecc0/6).*sin(3*f+2*aop0)...
    +3/8*s2i.*sin(4*f+2*aop0) + ecc0/16.*s2i.*sin(5*f+2*aop0));

dManS = A2./p.^2.*eta./ecc0.*(-(1-1.5*s2i).*((1-0.25*ecc0.^2).*sin(f) ...
    +0.5*ecc0.*sin(2*f) +ecc^2/12.*sin(3*f))...
    +s2i.*((0.25+5/16*ecc0.^2).*sin(f+2*aop0) - ecc0.^2/16.*sin(f-2*aop0)...
    -7/12.*(1-ecc0.^2/28).*sin(3*f+2*aop0) -3/8*ecc0.*sin(4*f+2*aop0)...
    -ecc0.^2/16.*sin(5*f+2*aop0)));
    
% Mean Short period variations    
c2f = (1-3*eta.^2+2*eta.^3)./ecc0.^2;

dEccM = A2./p.^2.*sin(inc0).^2.*(1-ecc0.^2)/6./ecc0.*c2f.*cos(2*aop0);

dIncM = -A2./p.^2/12.*sin(2*inc0).*c2f.*cos(2*aop0);

dRanM = -A2./p.^2/6.*cos(inc0).*c2f.*sin(2*aop0);

dAopM = A2./p.^2.*(s2i.*(1/8+(1-ecc0.^2)./ecc0.^2/6.*c2f) + c2f.*(1-s2i)/6).*sin(2*aop0);

dManM = -A2./p.^2.*eta.*s2i.*(1.8+(1/6./ecc0.^2+1/12).*c2f).*sin(2*aop0);

% Assign Elements
sma1 = sma0 + dSmaS;
ecc1 = ecc0 + dEccS - dEccM;
inc1 = inc0 + dIncS - dIncM;
ran1 = ran0 + dRanS - dRanM;
aop1 = aop0 + dAopS - dAopM;
man1 = man0 + dManS - dManM;

oeM = [sma1;ecc1;inc1;ran1;aop1;man1];