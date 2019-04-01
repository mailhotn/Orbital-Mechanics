% Orbits
inc = 24.6415;
ecc = 0;
nRepeats = 14;
nDays = 1;
sma = CalcRgtSma(ecc,inc,nRepeats,nDays);
p = sma*(1-ecc)^2;
primary = earth();
n = sqrt(primary.mu/sma^3);
eta = sqrt(1-ecc^2);

latEm = 10;
lonEm = 300;

% Ground Station Coordinates
meanAGs = wrapTo360([asind(sind(latEm)/sind(inc));
                     180 - asind(sind(latEm)/sind(inc))]);
raanGs = wrapTo360([lonEm - acosd(cosd(meanAGs(1))/cosd(latEm));
                    lonEm - acosd(cosd(meanAGs(2))/cosd(latEm))]);
                


% Create ROI
elevMin = 5;
% roiRad = 1;
roiRad = 90 - elevMin - asind(primary.Re/sma*sind(elevMin+90));
nPoints = 100;
z = cos(2*pi/nPoints) + 1i*sin(2*pi/nPoints);

meanRoi1 = nan(1,nPoints);
raanRoi1 = nan(1,nPoints);
meanRoi2 = nan(1,nPoints);
raanRoi2 = nan(1,nPoints);
for iPoint = 1:nPoints
    latPoint = sign(latEm + roiRad*imag(z^iPoint))*...
        min([abs(latEm + roiRad*imag(z^iPoint)),inc]);
    lonPoint = lonEm + roiRad*real(z^iPoint)/cosd(latPoint);
    
    meanRoi1(iPoint) = wrapTo360(asind(sind(latPoint)/sind(inc)));
    raanRoi1(iPoint) = wrapTo360(360+(lonPoint - ...
        acosd(cosd(meanRoi1(iPoint))/cosd(latPoint))));
    
    meanRoi2(iPoint) = wrapTo360(180-asind(sind(latPoint)/sind(inc)));
    raanRoi2(iPoint) = wrapTo360(360+(lonPoint - ...
        acosd(cosd(meanRoi2(iPoint))/cosd(latPoint))));
end
Area = polyarea(raanRoi1,meanRoi1);
% Create Constellation
latMat = [10,0;
          4,5];
Arch.nPlanes = latMat(1,1);
Arch.nAops = 1;
Arch.nSatsPerAop = latMat(2,2);
Phase.nC1 = latMat(2,1);
Phase.nC2 = 0;
Phase.nC3 = 0;
Orbit.inc = inc;
Orbit.ecc = 0;
Orbit.sma = sma;

InitCon = InitConDistOptimal(Arch,Phase,Orbit,latEm,lonEm);
% InitCon = InitConElliptical(ecc, inc, sma, latEm, lonEm);
raan0 = InitCon.raan1;
meanA0 = InitCon.M1;
[v,l] = eig(latMat);
nRelTraj = gcd(det(latMat),gcd(latMat(1,1)*nDays - latMat(1,2)*nRepeats,...
                               latMat(2,1)*nDays - latMat(2,2)*nRepeats));
nSats = det(latMat);
nSatsOrb = abs(gcd(latMat(1,2),latMat(2,2)));
nOrb = nSats/nSatsOrb;
meanA = nan(1,nSats);
raan = nan(1,nSats);
for iOrb = 1:nOrb
    for iSat = 1:nSatsOrb
        sol = latMat\[iOrb-1;iSat-1]*360;
        raan((iOrb-1)*nSatsOrb + iSat) = wrapTo360(sol(1));
        meanA((iOrb-1)*nSatsOrb + iSat) = wrapTo360(sol(2));
    end
end
raan = wrapTo360(raan + raan0);
meanA = wrapTo360(meanA + meanA0);
meanRate = 180/pi*(3/4*primary.J2*(primary.Re/p)^2*n*(5*cosd(inc)^2-1) +...
           3/4*primary.J2*(primary.Re/p)^2*n*eta*(3*cosd(inc)^2-1) + n);
raanRate = 180/pi*(-3/2*primary.J2*(primary.Re/p)^2*n*cosd(inc));
gsRate = primary.we;

rateVec = [raanRate - gsRate;
           meanRate];

figure(1)
clf
% Draw ROI
patch(raanRoi1,meanRoi1,[0.6,0.6,1])
hold on
roi = patch(raanRoi2,meanRoi2,[0.6,0.6,1]);

% Draw GS & Satellites
p = plot(raan,meanA,'o',raanGs,meanAGs,'r*','linewidth',2);
axis equal
xlim([0,360])
ylim([0,360])
xlabel('\Omega')
ylabel('\theta*')

% Draw Direction of movement
quiver(raan,meanA,rateVec(1)*ones(1,nSats),rateVec(2)*ones(1,nSats))
% quiver(raan,meanA,repmat(v(1,2)*100,1,nSats),repmat(v(2,2)*100,1,nSats),'linewidth',2)
legend([roi,p(1),p(2)],'ROI','Satellites','Ground Station')
hold off