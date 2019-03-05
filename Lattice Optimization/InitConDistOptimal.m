function InitCon = InitConDistOptimal(Arch,Phase,Orbit,latEm,lonEm)
% Create Constellation
nC = mod((Phase.nC1*Arch.nAops + Phase.nC3*(Arch.nSatsPerAop - Phase.nC2))/...
    gcd(Arch.nAops,Arch.nSatsPerAop - Phase.nC2),Arch.nPlanes);
nSatsOrb = (Arch.nAops*Arch.nSatsPerAop)/...
    gcd(Arch.nAops,Arch.nSatsPerAop - Phase.nC2);

latMat = [Arch.nPlanes 0;
    nC nSatsOrb];
nSats = det(latMat);
meanA = nan(1,nSats);
raan = nan(1,nSats);
for iOrb = 1:Arch.nPlanes
    for iSat = 1:nSatsOrb
        sol = latMat\[iOrb-1;iSat-1]*360;
        raan((iOrb-1)*nSatsOrb + iSat) = wrapTo360(sol(1));
        meanA((iOrb-1)*nSatsOrb + iSat) = wrapTo360(sol(2));
    end
end

% Find Emitter Location
meanAEm = wrapTo360([asind(sind(latEm)/sind(Orbit.inc));
    180 - asind(sind(latEm)/sind(Orbit.inc))]);
raanEm = wrapTo360([lonEm - acosd(cosd(meanAEm(1))/cosd(latEm));
    lonEm - acosd(cosd(meanAEm(2))/cosd(latEm))]);

% Optimize Initial Condition
ub = [360/nSatsOrb,360/Arch.nPlanes].';
options = optimoptions('fmincon','Display','off'); % For debugging
x0 = fmincon(@(x) DistCost(x,meanA,raan,meanAEm,raanEm),ub./2,eye(2),...
    360*ones(2,1),[],[],[0,0].',ub,[],options);
InitCon.M1 = x0(1);
InitCon.raan1 = x0(2);
InitCon.aop1 = 0;
end

function totDist = DistCost(x,meanA,raan,meanAEm,raanEm)
relPos1 = wrapTo360([meanA;raan]+x) - [meanAEm(1);raanEm(1)];
relPos2 = wrapTo360([meanA;raan]+x) - [meanAEm(2);raanEm(2)];
totDist = sum(min([dot(relPos1,relPos1,1),dot(relPos2,relPos2,1)]));
end