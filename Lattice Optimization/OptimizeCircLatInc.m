function [ Orbit ] = OptimizeCircLatInc(nRepeats,nDays,latGs,elevMin)
%OptimizeCircLatInc Finds optimal inclination from ROI area
%   Circular Orbits only
primary = earth();
inc = fminbnd(@(x) CalcRoiArea(x,nRepeats,nDays,latGs,elevMin,primary),latGs,89);
sma = CalcRgtSma(0,inc,nRepeats,nDays);
Orbit.ecc = 0;
Orbit.inc = inc;
Orbit.sma = sma;
end

function cost = CalcRoiArea(x,nRepeats,nDays,latGs,elevMin,primary)
sma = CalcRgtSma(0,x,nRepeats,nDays);
ecc = 0;
roiRad = 90 - elevMin - asind(primary.Re/sma*sind(elevMin+90));
nPoints = 100;
z = cos(2*pi/nPoints) + 1i*sin(2*pi/nPoints);

meanRoi1 = nan(1,nPoints);
raanRoi1 = nan(1,nPoints);
for iPoint = 1:nPoints
    latPoint = sign(latGs + roiRad*imag(z^iPoint))*...
        min([abs(latGs + roiRad*imag(z^iPoint)),x]);
    lonPoint = roiRad*real(z^iPoint)/cosd(latPoint);
    
    meanRoi1(iPoint) = (asind(sind(latPoint)/sind(x)));
    raanRoi1(iPoint) = (360 + (lonPoint - ...
        acosd(cosd(meanRoi1(iPoint))/cosd(latPoint))));
end
% %% Scaling by Rate
% p = sma*(1-ecc)^2;
% eta = sqrt(1-ecc^2);
% n = sqrt(primary.mu/sma^3);
% meanRate = 180/pi*(3/4*primary.J2*(primary.Re/p)^2*n*(5*cosd(x)^2-1) +...
%            3/4*primary.J2*(primary.Re/p)^2*n*eta*(3*cosd(x)^2-1) + n);
% raanRate = 180/pi*(-3/2*primary.J2*(primary.Re/p)^2*n*cosd(x));
% gsRate = primary.we;
% rateVec = [raanRate - gsRate;
%            meanRate];
% scale = abs(rateVec(1)/rateVec(2));
% meanRoi1 = meanRoi1*scale;
%% Scaled Area
cost = -polyarea(raanRoi1,meanRoi1);
end