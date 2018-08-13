function [ InitCon , t2 ] = InitConElliptical( ecc, inc, sma, latGs, lonGs, gmst0 )
%InitConTdoaElliptical returns the initial aop, raan & mean anomaly
% for a satellite in an elliptical orbit such that on the first orbit the
% satellite passes above a ground station at apogee
if nargin < 6
    gmst0 = 0;
end
InitCon = struct();
primary = earth;

aol1 = 0;
t1 = 0;
% Rates are in rad/s, final results in deg
n = sqrt(primary.mu/sma^3); % rad/s
p = sma*(1-ecc^2);
J2Mean   = 3/4*primary.J2*(primary.Re/p)^2; % rad/s
aopRate  = -J2Mean*n*(1 - 5*cosd(inc)^2);
meanRate = n*(1-J2Mean*sqrt(1-ecc^2)*(1 - 3*cosd(inc)^2));
raanRate = -2*J2Mean*n*cosd(inc);

aol2 = asin(sind(latGs)/sind(inc));
t2 = (aol2 - aol1)/(meanRate + aopRate) + t1;

InitCon.aop1 = 180/pi*(aol2 - pi - aopRate*(t2-t1));
InitCon.M1 = 180/pi*(pi - meanRate*(t2-t1));
InitCon.raan1 = gmst0 + lonGs + (primary.we - 180/pi*raanRate)*(t2-t1) - ...
    atan2d(cosd(inc),(tan(aol2-aol1))^-1);
end
