function [ Orbit ] = OptimizeCircLatInc(nRepeats,nDays,latGs,elevMin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
primary = earth();
inc = fminbnd(@(x) CalcRoiArea(x,nRepeats,nDays,latGs,elevMin,primary),latGs,89);
sma = CalcRgtSma(0,inc,nRepeats,nDays);
Orbit.ecc = 0;
Orbit.inc = inc;
Orbit.sma = sma;
end

function cost = CalcRoiArea(x,nRepeats,nDays,latGs,elevMin,primary)
sma = CalcRgtSma(0,x,nRepeats,nDays);
roiRad = 90 - elevMin - asind(primary.Re/sma*sind(elevMin+90));
nPoints = 100;
z = cos(2*pi/nPoints) + 1i*sin(2*pi/nPoints);

meanRoi1 = nan(1,nPoints);
raanRoi1 = nan(1,nPoints);
for iPoint = 1:nPoints
    latPoint = min([latGs + roiRad*imag(z^iPoint),x]);
    lonPoint = roiRad*real(z^iPoint);
    
    meanRoi1(iPoint) = (asind(sind(latPoint)/sind(x)));
    raanRoi1(iPoint) = (360 + (lonPoint - ...
        acosd(cosd(meanRoi1(iPoint))/cosd(latPoint))));
end
cost = -polyarea(raanRoi1,meanRoi1);
end