%% Define Scenarios for Analysis
clear
paretoFolder = 'C:\Users\User\Google Drive\Master''s Degree\Lattice Optimization Data\Final Results\OG Pareto Lists';


latList = [30,40,50,60];
%% Walker
ConFolder = 'C:\Users\User\Google Drive\Master''s Degree\Walker Optimization Data\Previous Optimization Runs\Latticified Walker multilat multiInc collision';

ConTable = [];
for iLat = 1:length(latList)
    load([paretoFolder '\ParetoList_Lat_' num2str(latList(iLat)) '_Walker Ex.mat'])
    for iPar = 1:length(paretoSet.nSats)
        load([ConFolder '\LatticeExSol_Lat_' num2str(latList(iLat))...
            '_nSats_' num2str(paretoSet.nSats(iPar)) '_hA_' ...
            num2str(paretoSet.hAs(iPar)) '.mat']);
        iCon = find(ExSol.archMat(1,:) == paretoSet.nPlanes(iPar));
        Orbit = ExSol.orbits{iCon};
        nGroundTracks = ExSol.archMat(2,iCon)*ExSol.archMat(1,iCon)...
            /gcd(ExSol.archMat(1,iCon),ExSol.phaseMat(3,iCon));
        paramVec = [latList(iLat) ExSol.nSats ExSol.archMat(1,iCon).' -ExSol.phaseMat(1,iCon).'...
            Orbit.sma Orbit.inc nGroundTracks];
        ConTable = [ConTable; paramVec];
    end
end
matrix2lyx(ConTable,'WalkerParameters.lyx')


%% LFC Exhaustive
ConFolder = 'C:\Users\User\Google Drive\Master''s Degree\Lattice Optimization Data\Previous Runs\LatticeDef v1';
ConTable = [];
% conVec1 = [];
% conVec2 = [];
for iLat = 1:length(latList)
    load([paretoFolder '\ParetoList_Lat_' num2str(latList(iLat)) '_LFC Ex.mat'])
    for iPar = 1:length(paretoSet.nSats)
        load([ConFolder '\LatticeExSol_Lat_' num2str(latList(iLat))...
            '_nSats_' num2str(paretoSet.nSats(iPar)) '_hA_' ...
            num2str(paretoSet.hAs(iPar)) '.mat']);
        iCon = find(ExSol.archMat(1,:) == paretoSet.nPlanes(iPar));
        Orbit = ExSol.orbits{iCon};
        %         conVec1 = [conVec1; ~mod(ExSol.archMat(1,iCon),ExSol.phaseMat(3,iCon))];
        %         conVec2 = [conVec2; ~(((14-ExSol.phaseMat(1,iCon))/...
        %         gcd(ExSol.archMat(1,iCon),ExSol.phaseMat(3,iCon)))-...
        %         round((14-ExSol.phaseMat(1,iCon))/...
        %         gcd(ExSol.archMat(1,iCon),ExSol.phaseMat(3,iCon))))];
        nGroundTracks = ExSol.archMat(2,iCon)*ExSol.archMat(1,iCon)...
            /gcd(ExSol.archMat(1,iCon),ExSol.phaseMat(3,iCon));
        paramVec = [latList(iLat) ExSol.nSats ExSol.archMat(:,iCon).' ExSol.phaseMat(:,iCon).'...
            Orbit.sma Orbit.inc Orbit.hA nGroundTracks];
        ConTable = [ConTable; paramVec];
    end
end
matrix2lyx(ConTable,'LatticeExParameters.lyx')

%% GA
ConFolder = 'C:\Users\User\Google Drive\Master''s Degree\Lattice Optimization Data\GA Standard\Previous Runs\Version 6 - definitive';
ConTable = [];
conVec1 = []; % Nc3 constraint
conVec2 = []; % Nc1 constraint
for iLat = 1:length(latList)
    load([paretoFolder '\ParetoList_Lat_' num2str(latList(iLat)) '_LFC Ga.mat'])
    for iPar = 1:length(paretoSet.nSats)
        load([ConFolder '\LatticeGaSol_Lat_' num2str(latList(iLat))...
            '_nSats_' num2str(paretoSet.nSats(iPar)) '.mat']);
        iCon = find(GaSol.archMat(1,:) == paretoSet.nPlanes(iPar));
        Orbit = GaSol.orbits{iCon};
        conVec1 = [conVec1; ~mod(GaSol.archMat(1,iCon),GaSol.phaseMat(3,iCon))];
%         conVec1 = [conVec1; (gcd(GaSol.archMat(1,iCon),GaSol.phaseMat(3,iCon))~=1) || GaSol.phaseMat(3,iCon)==1];
        conVec2 = [conVec2; ~(((14-GaSol.phaseMat(1,iCon))/...
        gcd(GaSol.archMat(1,iCon),GaSol.phaseMat(3,iCon)))-...
        round((14-GaSol.phaseMat(1,iCon))/...
        gcd(GaSol.archMat(1,iCon),GaSol.phaseMat(3,iCon))))];
          
        paramVec = [latList(iLat) GaSol.nSats GaSol.archMat(:,iCon).' GaSol.phaseMat(:,iCon).'...
            Orbit.sma Orbit.inc Orbit.hA];
        ConTable = [ConTable; paramVec];
    end
end
matrix2lyx([ConTable conVec1.*conVec2],'LatticeGaParameters.lyx')

%% LFC Exhaustive by hA
ConFolder = 'C:\Users\User\Google Drive\Master''s Degree\Lattice Optimization Data\Previous Runs\LatticeDef v1';
ConTable = [];
iHa = 2;
haList = [0, 900, 1000];
% conVec1 = [];
% conVec2 = [];
for iLat = 1:length(latList)
    load([paretoFolder '\ParetoList_Lat_' num2str(latList(iLat)) '_LFC Ex_hASplit.mat'])
    for iPar = 1:length(paretoSet.nSats{iHa})
        load([ConFolder '\LatticeExSol_Lat_' num2str(latList(iLat))...
            '_nSats_' num2str(paretoSet.nSats{iHa}(iPar)) '_hA_' ...
            num2str(haList(iHa)) '.mat']);
        iCon = find(ExSol.archMat(1,:) == paretoSet.nPlanes{iHa}(iPar));
        Orbit = ExSol.orbits{iCon};
    
        nGroundTracks = ExSol.archMat(2,iCon)*ExSol.archMat(1,iCon)...
            /gcd(ExSol.archMat(1,iCon),ExSol.phaseMat(3,iCon));
        paramVec = [latList(iLat) ExSol.nSats ExSol.archMat(:,iCon).' ExSol.phaseMat(:,iCon).'...
            Orbit.sma Orbit.inc Orbit.hA nGroundTracks];
        ConTable = [ConTable; paramVec];
    end
end
matrix2lyx(ConTable,'LatticeExParameters.lyx')
