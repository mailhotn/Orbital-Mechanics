function [ Orbit ] = CalcOptIncRoiA(nRepeats,nDays,latEm,hA,elevMin)
%CalcOptIncRoiA Finds optimal inclination from ROI area
%   Outputs Orbit struct for elliptical constellation
primary = earth();
inc = fminbnd(@(x) CalcRoiArea(x,nRepeats,nDays,latEm,elevMin,primary),latEm,89);
[sma,ecc] = CalcRgtSmaApoHeight(inc,hA,nRepeats,nDays);
Orbit.ecc = ecc;
Orbit.inc = inc;
Orbit.sma = sma;
Orbit.hA = hA;
end

function cost = CalcRoiArea(x,nRepeats,nDays,latEm,elevMin,primary)
sma = CalcRgtSma(0,x,nRepeats,nDays);
roiRad = 90 - elevMin - asind(primary.Re/sma*sind(elevMin+90));
nPoints = 100;
z = cos(2*pi/nPoints) + 1i*sin(2*pi/nPoints);

meanRoi1 = nan(1,nPoints);
raanRoi1 = nan(1,nPoints);
for iPoint = 1:nPoints
    %     latPoint = sign(latEm + roiRad*imag(z^iPoint))*...
    %         min([abs(latEm + roiRad*imag(z^iPoint)),x]);
    %     lonPoint = roiRad*real(z^iPoint)/cosd(latPoint);
    
    phiE = wrapTo360(angle(z^iPoint)*180/pi);
    latPEm = 90 - latEm;
    latPP = acosd(cosd(roiRad)*cosd(latPEm) +...
        sind(roiRad)*sind(latPEm)*cosd(phiE));
    latPoint = 90 - latPP;
    latPoint = sign(latPoint)*min([abs(latPoint),x]);
    latPP = 90 - latPoint;
    dLon = real(acosd((cosd(roiRad) - cosd(latPP)*cosd(latPEm))/...
        (hemid(latPEm)*sind(latPEm)*sind(latPP))) - 90*(hemid(latPEm)-1));
    lonPoint = dLon*hemid(phiE);
    
    meanRoi1(iPoint) = (asind(sind(latPoint)/sind(x)));
    raanRoi1(iPoint) = (360 + (lonPoint - ...
        acosd(cosd(meanRoi1(iPoint))/cosd(latPoint))));
end
% %% Scaling by Rate
% ecc = 0;
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
% meanRoi1 = meanRoi1/rateVec(2);
% raanRo11 = raanRoi1/abs(rateVec(1));
%%
cost = -polyarea(raanRoi1,meanRoi1);
end