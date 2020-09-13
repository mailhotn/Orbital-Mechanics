clear
dataFolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Final Results';
nameList = {'Walker Ex','LFC Ex','LFC Ga'};
folderList = {...
    'C:\Users\User\Dropbox\Walker Optimization Data\Previous Optimization Runs\Latticified Walker multilat multiInc collision';...
    'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\LatticeDef v1';...
    'C:\Users\User\Dropbox\Lattice Optimization Data\GA Standard\Previous Runs\Version 6 - definitive';...
    };
latEm = 30;
hAList = [0,900,1000];

WalkerExPhase = [];
LatticeExPhase = cell(3,1);
LatticeGaPhase = [];

WalkerExOrbit = [];
LatticeExOrbit = cell(3,1);
LatticeGaOrbit = [];

WalkerExArch = [];
LatticeExArch = cell(3,1);
LatticeGaArch = [];

WalkerMinDist = [];
LatticeExMinDist = cell(3,1);
LatticeGaMinDist = [];

for iScenario = 2
    
    load([dataFolder '/ParetoList_Lat_' num2str(latEm) '_' nameList{iScenario} '_hASplit.mat']);
    if iScenario == 2
        hAList = [0, 900, 1000];
    else
        hAList = 0;
    end
    for iHA = 3:length(hAList)
        for iCon = 1:length(paretoSet.nPlanes{iHA})
            switch iScenario
                case 1 % Walker
                    load([folderList{iScenario} '\LatticeExSol_Lat_'...
                        num2str(latEm) '_nSats_' num2str(paretoSet.nSats{iHA}(iCon)) ...
                        '_hA_' num2str(hAList(iHA)) '.mat']);
                    iPar = find(ExSol.archMat(1,:) == paretoSet.nPlanes{iHA}(iCon));
                    
                    WalkerExPhase = [WalkerExPhase, ExSol.phaseMat(1:3,iPar)];
                    WalkerExOrbit = [WalkerExOrbit, [ExSol.orbits{iPar}.sma;...
                        ExSol.orbits{iPar}.ecc;ExSol.orbits{iPar}.inc;]];
                    WalkerExArch = [WalkerExArch, [paretoSet.nSats{iHA}(iCon);...
                        paretoSet.nPlanes{iHA}(iCon)]];
                    
                    Arch = struct();
                    Arch.nSats = paretoSet.nSats{iHA}(iCon);
                    Arch.nPlanes = ExSol.archMat(1,iPar);
                    Arch.nAops = ExSol.archMat(2,iPar);
                    Arch.nSatsPerAop = ExSol.archMat(3,iPar);
                    
                    Phase = struct();
                    Phase.nC1 = ExSol.phaseMat(1,iPar);
                    Phase.nC2 = ExSol.phaseMat(2,iPar);
                    Phase.nC3 = ExSol.phaseMat(3,iPar);
                    
                    Orbit = ExSol.orbits{iPar};
                    
                    InitCon = ExSol.inits{iPar};
                    
                case 2 % LFC Ex
                    load([folderList{iScenario} '\LatticeExSol_Lat_'...
                        num2str(latEm) '_nSats_' num2str(paretoSet.nSats{iHA}(iCon)) ...
                        '_hA_' num2str(hAList(iHA)) '.mat']);
                    
                    iPar = find(ExSol.archMat(1,:) == paretoSet.nPlanes{iHA}(iCon));
                    
                    LatticeExPhase{iHA} = [LatticeExPhase{iHA}, ExSol.phaseMat(1:3,iPar)];
                    LatticeExOrbit{iHA} = [LatticeExOrbit{iHA}, [ExSol.orbits{iPar}.sma;...
                        ExSol.orbits{iPar}.ecc;ExSol.orbits{iPar}.inc;]];
                    LatticeExArch{iHA} = [LatticeExArch{iHA}, [paretoSet.nSats{iHA}(iCon);...
                        paretoSet.nPlanes{iHA}(iCon)]];
                    
                    Arch = struct();
                    Arch.nSats = paretoSet.nSats{iHA}(iCon);
                    Arch.nPlanes = ExSol.archMat(1,iPar);
                    Arch.nAops = ExSol.archMat(2,iPar);
                    Arch.nSatsPerAop = ExSol.archMat(3,iPar);
                    
                    Phase = struct();
                    Phase.nC1 = ExSol.phaseMat(1,iPar);
                    Phase.nC2 = ExSol.phaseMat(2,iPar);
                    Phase.nC3 = ExSol.phaseMat(3,iPar);
                    
                    Orbit = ExSol.orbits{iPar};
                    
                    InitCon = ExSol.inits{iPar};
                    
                case 3 % LFC GA
                    load([folderList{iScenario} '\LatticeGaSol_Lat_' num2str(latEm)...
                        '_nSats_' num2str(paretoSet.nSats{iHA}(iCon)) '.mat']);
                    iPar = find(GaSol.archMat(1,:) == paretoSet.nPlanes{iHA}(iCon));
                    
                    LatticeGaPhase = [LatticeGaPhase, GaSol.phaseMat(1:3,iPar)];
                    LatticeGaOrbit = [LatticeGaOrbit, [GaSol.orbits{iPar}.sma;...
                        GaSol.orbits{iPar}.ecc;GaSol.orbits{iPar}.inc;]];
                    LatticeGaArch = [LatticeGaArch, [paretoSet.nSats{iHA}(iCon);...
                        paretoSet.nPlanes{iHA}(iCon)]];
                    
                    Arch = struct();
                    Arch.nSats = paretoSet.nSats{iHA}(iCon);
                    Arch.nPlanes = GaSol.archMat(1,iPar);
                    Arch.nAops = GaSol.archMat(2,iPar);
                    Arch.nSatsPerAop = GaSol.archMat(3,iPar);
                    
                    Phase = struct();
                    Phase.nC1 = GaSol.phaseMat(1,iPar);
                    Phase.nC2 = GaSol.phaseMat(2,iPar);
                    Phase.nC3 = GaSol.phaseMat(3,iPar);
                    
                    Orbit = GaSol.orbits{iPar};
                    
                    InitCon = GaSol.inits{iPar};
                    
            end
            
            
            LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);
            minDist = CalcMinDist(LC);
%             minDist = minDist/LC.sma;
            switch iScenario
                case 1
                    WalkerMinDist = [WalkerMinDist, minDist];
                case 2
                    LatticeExMinDist{iHA} = [LatticeExMinDist{iHA}, minDist];
                case 3
                    LatticeGaMinDist = [LatticeGaMinDist, minDist];
            end
            
                    figure(iScenario*100 + iCon*10 + 1)
            %         subplot(2,1,1)
            %         PlotGroundTrack(LC);
            %         title([nameList{iScenario} ' ' num2str(Arch.nSats) '/' ...
            %             num2str(Arch.nPlanes) ' ' num2str(Orbit.hA)])
            %
            %         subplot(2,1,2)
                    PlotPdopMap(LC,latEm,0,5,true,true);
        end
    end
end

%% Generate example of too regional coverage
clear
close all
dataFolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\LatticeDef v1';
load([dataFolder '\LatticeExSol_Lat_30_nSats_60_hA_900.mat']);

Arch.nSats = 60;
Arch.nPlanes = 60;
Arch.nAops = 1;
Arch.nSatsPerAop = 1;

Orbit = ExSol.orbits{end};

Phase.nC1 = ExSol.phaseMat(1,end);
Phase.nC2 = ExSol.phaseMat(2,end);
Phase.nC3 = ExSol.phaseMat(3,end);

InitCon = ExSol.inits{end};

LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);
time = 0:10:86164;
Prop = Propagator(LC);
[~, xEci] = Prop.PropEciJ2(time);

primary = earth();
gmst = primary.we*time;
nSats = size(xEci,2)/6;
tracks = nan(2*nSats,length(gmst));
for iTime = 1:length(gmst)
    xEcef = eci2ecef(reshape(xEci(iTime,:).',6,nSats),gmst(iTime));
    rEcef = xEcef(1:3,:)./sqrt(dot(xEcef(1:3,:),xEcef(1:3,:),1));
    lla   = ecef2lla(rEcef.'*primary.Re,0,primary.Re).';
    tracks(:,iTime) = reshape(lla(1:2,:),2*nSats,1);
end

figure(1)
PlotPdopMap(LC,30,0,5,true);
hold on
for iSat = 1:nSats
        plot(wrapTo180(tracks(2*iSat,:)),tracks(2*iSat-1,:),'.w');
end
hold off
    

