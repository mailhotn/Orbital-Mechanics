% Orbits
inc = 50;
ecc = 0;
nRepeats = 14;
nDays = 1;
sma = CalcRgtSma(ecc,inc,nRepeats,nDays);
p = sma*(1-ecc)^2;
primary = earth();
n = sqrt(primary.mu/sma^3);
eta = sqrt(1-ecc^2);

latGs = 40;
lonGs = 0;

dLat = 10;
dLon = 10;

meanAGs = wrapTo360([asind(sind(latGs)/sind(inc));
                     180 - asind(sind(latGs)/sind(inc))]);
raanGs = wrapTo360([lonGs - acosd(cosd(meanAGs(1))/cosd(latGs));
                    lonGs - acosd(cosd(meanAGs(2))/cosd(latGs))]);
meanRoi1 = nan(1,8);
raanRoi1 = nan(1,8);
meanRoi2 = nan(1,8);
raanRoi2 = nan(1,8);
dLatList = min([-1,-1,-1,0,1,1,1,0]*dLat,ones(1,8)*(inc-latGs));
dLonList = [-1,0,1,1,1,0,-1,-1]*dLon;
for iPoint = 1:8
    meanRoi1(iPoint) = (asind(sind(latGs + dLatList(iPoint))/sind(inc)));
    raanRoi1(iPoint) = 360+(lonGs + dLonList(iPoint) - ...
        acosd(cosd(meanRoi1(iPoint))/cosd(latGs + dLatList(iPoint))));
    
    meanRoi2(iPoint) = (180-asind(sind(latGs + dLatList(iPoint))/sind(inc)));
    raanRoi2(iPoint) = 360+(lonGs + dLonList(iPoint) - ...
        acosd(cosd(meanRoi2(iPoint))/cosd(latGs + dLatList(iPoint))));
end

latMat = [13,0;
          6,5];
raan0 = -40.6813;
meanA0 = 0;
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
% quiver(0,0,v(1,2)*100,v(2,2)*100)
legend([roi,p(1),p(2)],'ROI','Satellites','Ground Station')
hold off