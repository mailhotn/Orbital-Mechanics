% Given a folder of Walker RGT results, this script saves all results as
% equivalent Lattice Flower Constellations
%% Define Folder to reformat
walkerFolder = ['C:\Users\User\Dropbox\Walker Optimization Data'...
    '\Previous Optimization Runs\Walker RGT Ex Search delta inc 10'];
latList = 10:10:80;
minSats = 20;
maxSats = 80;
nSats = minSats:maxSats;
nCons = numel(nSats);
primary = earth();
%% Reformat Files as Lattice Flower Constellations
for iLat = 1:numel(latList)
    for iSats = 1:nCons
        load([walkerFolder '\WalkerRgtExSol_Lat_' num2str(latList(iLat))...
            '_T_' num2str(nSats(iSats)) '.mat']);
        WSol = ExSol;
        clear ExSol;
        
        [intVec,indOpt] = min(WSol.intPdop(:,divisors(WSol.nSats)));
        maxVec = WSol.maxPdop(sub2ind(size(WSol.maxPdop),indOpt,...
            divisors(WSol.nSats)));
        covVec = WSol.coverage(sub2ind(size(WSol.coverage),indOpt,...
            divisors(WSol.nSats)));
        [fit,iOpt] = min(intVec);
        
        phaseMat = [-(indOpt-1);zeros(2,length(indOpt))];
        archMat = [divisors(WSol.nSats);
            ones(1,length(divisors(WSol.nSats)));
            WSol.nSats./divisors(WSol.nSats)];
        
        Orbit.ecc = 0;
        Orbit.sma = primary.Re + WSol.alt;
        Orbit.inc = WSol.inc;
        
        InitCon.M1 = 0;
        InitCon.raan1 = WSol.raan0;
        InitCon.aop1 = 0;
        
        ExSol.nSats = WSol.nSats;
        ExSol.archMat = archMat;
        ExSol.phaseMat = phaseMat;
        ExSol.Orbit = Orbit;
        ExSol.InitCon = InitCon;
        ExSol.latGs = WSol.latGs;
        ExSol.optNPlanes = archMat(1,iOpt);
        ExSol.iOpt = iOpt;
        ExSol.PropParams = WSol.PropParams;
        ExSol.coverage = covVec;
        ExSol.maxPdop  = maxVec;
        ExSol.intPdop  = intVec;
        ExSol.fit = fit;
        
        save([walkerFolder '\LatticeExSol_Lat_' num2str(latList(iLat))...
            '_nSats_' num2str(nSats(iSats)) '_ecc_0.mat'], 'ExSol');
        clear ExSol;
    end
end

%% Creat OptParams Struct
clear
PropParams.maxPdop = 1000;
PropParams.timeVec = 0:10:86164;
PropParams.elevMin = 5;
PropParams.relTol = 1e-6;
PropParams.absTol = 1e-6;
PropParams.datafolder = ['C:\Users\User\Dropbox\Walker Optimization Data'...
    '\Previous Optimization Runs\Walker RGT Ex Search delta inc 10'];

primary = earth();
lonGs = 0;

nRepeats = 14;
nDays = 1;

latList = 10:10:80;
maxSats = 80;
minSats = 20;

delInc = 10;
eccList = 0;

save([PropParams.datafolder '\OptParams.mat'])