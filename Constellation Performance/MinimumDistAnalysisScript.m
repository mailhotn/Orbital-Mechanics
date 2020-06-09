clear
datafolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\Apogee Height x3, del inc 10, multilat PDOP';
load([datafolder '\OptParams.mat']);

lat = 60;
nSats = 60;
nPlanes = 6;
hA = 900;
load([datafolder '\LatticeExSol_Lat_' num2str(lat) '_nSats_' num2str(nSats)...
    '_hA_' num2str(hA) '.mat']);

%% Determine Minimum Distance

iCon = ExSol.archMat(1,:) == nPlanes;
Arch.nPlanes = ExSol.archMat(1,iCon);
Arch.nAops = ExSol.archMat(2,iCon);
Arch.nSatsPerAop = ExSol.archMat(3,iCon);

Phase.nC1 = ExSol.phaseMat(1,iCon);
Phase.nC2 = ExSol.phaseMat(2,iCon);
Phase.nC3 = ExSol.phaseMat(3,iCon);


LC = LatticeConstellation(Arch,Phase,ExSol.Orbit,ExSol.InitCon);
Con = LC;
prop = Propagator(Con);

timeVec = 0:10:86164;
[~, xEci] = prop.PropEciJ2(timeVec);
xEci = reshape(xEci.',6,Con.nSats*length(timeVec));
rEci = xEci(1:3,:);
rEci = reshape(rEci,3*Con.nSats,length(timeVec));
xEci = reshape(xEci,Con.nSats*6,length(timeVec));
xEci = xEci.';
minDist = inf;
for iSat = 1:Con.nSats
    satsList = iSat+1:Con.nSats;
    for iSat2 = satsList
        distVecs = rEci(3*iSat-2:3*iSat,:)-rEci(3*iSat2-2:3*iSat2,:);
        dist = sqrt(dot(distVecs,distVecs,1));
        if min(dist) < minDist
            [minDist,minTime] = min(dist);
            minSats = [iSat,iSat2];
        end
    end
end
