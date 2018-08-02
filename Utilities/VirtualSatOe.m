function [ inc, sma, raan0 ] = VirtualSatOe( ecc, Rgt, GroundS, nOrbits, primary )
%VirtualSatOe calculates the inc, sma and raan0 for an RGT constellation
%   Based on procedure in T. Shtark and P. Gurfil,
%  “Regional positioning using a low Earth orbit satellite constellation,” 2018
if nargin < 5
    primary = earth();
end
inc = fminbnd(@(x)VirtualSatCost(x, ecc, Rgt, GroundS, nOrbits, primary),0,180);

sma = CalcRgtSma(ecc, inc, Rgt.jRepeats, Rgt.kDays, primary);

n = sqrt(primary.mu/sma^3);
raanDot = -1.5*(primary.J2*sqrt(primary.mu)*primary.Re^2)/(sma^3.5*(1-ecc^2)^2)*cosd(inc);
aol0 = 0;
t0 = 0;
aol1 = asin(sind(GroundS.lat)/sind(inc));
t1 = (aol1-aol0)/n + t0;

raan0 = GroundS.lon + (primary.we - raanDot*180/pi)*(t1-t0) - ...
        atan2d(cosd(inc),(tan(aol1 - aol0))^-1);

end

function [ f ] = VirtualSatCost( x, ecc, Rgt, GroundS, nOrbits, primary )
%VirtualSatCost is the function to minimize when solving for a virtual sat
sma = CalcRgtSma(ecc, x, Rgt.jRepeats, Rgt.kDays, primary);
n = sqrt(primary.mu/sma^3);
raanDot = -1.5*(primary.J2*sqrt(primary.mu)*primary.Re^2)/(sma^3.5*(1-ecc^2)^2)*cosd(x);
aol0 = 0;
t0 = 0;
aol1 = asin(sind(GroundS.lat)/sind(x));
t1 = (aol1-aol0)/n + t0;
aol2 = pi - aol1;
t2 = (nOrbits*2*pi+aol2-aol0)/n + t0;

F = -(primary.we*pi/180 - raanDot)*(t2-t1) - ...
    atan2(cosd(x),(tan(aol1 - aol0))^-1) + ...
    atan2(cosd(x),(tan(aol2 - aol0))^-1);
f = F^2;
end

