%% Load Optimization Parameters

lat10Folder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\Lattice Version 3';
latIncFolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\Lattice O-M Optimal Inc + Start v1';
latOptFolder = 'C:\Users\User\Dropbox\Lattice Optimization Data';

walkerFolder = ['C:\Users\User\Dropbox\Walker Optimization Data'...
    '\Previous Optimization Runs\Walker RGT Ex Search delta inc 10'];

load ([latOptFolder '\OptParams.mat']);
% OptParams.mat contains the following variables and structures:
%
% PropParams.maxPdop
% PropParams.timeVec
% PropParams.elevMin
% PropParams.relTol
% PropParams.absTol
% PropParams.datafolder
%
% primary
% lonGs
%
% nRepeats
% nDays
%
% latList
% maxSats
% minSats
%
% delInc
% eccList

%% Initialize Data Matrices
nLats = length(latList);
nCons = maxSats - minSats + 1;
nEcc = length(eccList);

lat10Perf.maxPdop  = nan(nEcc,nCons,nLats);
lat10Perf.intPdop  = nan(nEcc,nCons,nLats);
lat10Perf.coverage = nan(nEcc,nCons,nLats);
latIncPerf.maxPdop  = nan(nEcc,nCons,nLats);
latIncPerf.intPdop  = nan(nEcc,nCons,nLats);
latIncPerf.coverage = nan(nEcc,nCons,nLats);
latOptPerf.maxPdop  = nan(nEcc,nCons,nLats);
latOptPerf.intPdop  = nan(nEcc,nCons,nLats);
latOptPerf.coverage = nan(nEcc,nCons,nLats);
walkerPerf.maxPdop   = nan(1,nCons,nLats);
walkerPerf.intPdop   = nan(1,nCons,nLats);
walkerPerf.coverage  = nan(1,nCons,nLats);
% inc       = nan(nLats,nCons);
% alt       = nan(nLats,nCons);
% latticeNSats    = repmat(minSats:maxSats,nLats,1);
nSats           = repmat(minSats:maxSats,nEcc,1);

%% Performance Goals
intTarget = 5;
maxTarget = 20;
covTarget = 99.5;
lat10Perf.nSatsToAchieve = nan(nEcc,nLats);
lat10Perf.nPlanesToAchieve = nan(nEcc,nLats);
latOptPerf.nSatsToAchieve = nan(nEcc,nLats);
latIncPerf.nPlanesToAchieve = nan(nEcc,nLats);
latIncPerf.nSatsToAchieve = nan(nEcc,nLats);
latOptPerf.nPlanesToAchieve = nan(nEcc,nLats);
walkerPerf.nSatsToAchieve = nan(1,nLats);
walkerPerf.nPlanesToAchieve = nan(1,nLats);

%% Check Performance
for iLat = 1:nLats
    %     latticeMaxPdop  = nan(nEcc,nCons);
    %     latticeIntPdop  = nan(nEcc,nCons);
    %     latticeCoverage = nan(nEcc,nCons);
    %     walkerMaxPdop   = nan(1,nCons);
    %     walkerIntPdop   = nan(1,nCons);
    %     walkerCoverage  = nan(1,nCons);
    for iEcc = 1:length(eccList)
        for iSats = 1:nCons
            
            % Lattice del inc = 10
            load([lat10Folder '\LatticeExSol_Lat_' num2str(latList(iLat))...
                '_nSats_' num2str(nSats(iEcc,iSats)) '_ecc_' num2str(eccList(iEcc)) '.mat']);
            lat10Perf.maxPdop(iEcc,iSats,iLat) = ExSol.maxPdop(ExSol.iOpt);
            lat10Perf.intPdop(iEcc,iSats,iLat) = ExSol.intPdop(ExSol.iOpt);
            lat10Perf.coverage(iEcc,iSats,iLat) = ExSol.coverage(ExSol.iOpt);
            % Check Goal
            if isnan(lat10Perf.nSatsToAchieve(iEcc,iLat))
                achieveInt = ExSol.intPdop < intTarget;
                achieveMax = ExSol.maxPdop < maxTarget;
                achieveCov = ExSol.coverage > covTarget;
                if any(achieveInt & achieveMax & achieveCov)
                    lat10Perf.nSatsToAchieve(iEcc, iLat) = nSats(iEcc,iSats);
                    planeVec = achieveInt & achieveMax & achieveCov;
                    [~,iMinPlanes] = min(~(planeVec>0));
                    lat10Perf.nPlanesToAchieve(iEcc,iLat) = ExSol.archMat(1,iMinPlanes);
                end
            end
            
            % Lattice Optimal inc
            load([latIncFolder '\LatticeExSol_Lat_' num2str(latList(iLat))...
                '_nSats_' num2str(nSats(iEcc,iSats)) '_ecc_' num2str(eccList(iEcc)) '.mat']);
            latIncPerf.maxPdop(iEcc,iSats,iLat) = ExSol.maxPdop(ExSol.iOpt);
            latIncPerf.intPdop(iEcc,iSats,iLat) = ExSol.intPdop(ExSol.iOpt);
            latIncPerf.coverage(iEcc,iSats,iLat) = ExSol.coverage(ExSol.iOpt);
            % Check Goal
            if isnan(latIncPerf.nSatsToAchieve(iEcc,iLat))
                achieveInt = ExSol.intPdop < intTarget;
                achieveMax = ExSol.maxPdop < maxTarget;
                achieveCov = ExSol.coverage > covTarget;
                if any(achieveInt & achieveMax & achieveCov)
                    latIncPerf.nSatsToAchieve(iEcc, iLat) = nSats(iEcc,iSats);
                    planeVec = achieveInt & achieveMax & achieveCov;
                    [~,iMinPlanes] = min(~(planeVec>0));
                    latIncPerf.nPlanesToAchieve(iEcc,iLat) = ExSol.archMat(1,iMinPlanes);
                end
            end
            
            % Lattice Optimal inc & Start
            load([latOptFolder '\LatticeExSol_Lat_' num2str(latList(iLat))...
                '_nSats_' num2str(nSats(iEcc,iSats)) '_ecc_' num2str(eccList(iEcc)) '.mat']);
            latOptPerf.maxPdop(iEcc,iSats,iLat) = ExSol.maxPdop(ExSol.iOpt);
            latOptPerf.intPdop(iEcc,iSats,iLat) = ExSol.intPdop(ExSol.iOpt);
            latOptPerf.coverage(iEcc,iSats,iLat) = ExSol.coverage(ExSol.iOpt);
            % Check Goal
            if isnan(latOptPerf.nSatsToAchieve(iEcc,iLat))
                achieveInt = ExSol.intPdop < intTarget;
                achieveMax = ExSol.maxPdop < maxTarget;
                achieveCov = ExSol.coverage > covTarget;
                if any(achieveInt & achieveMax & achieveCov)
                    latOptPerf.nSatsToAchieve(iEcc, iLat) = nSats(iEcc,iSats);
                    planeVec = achieveInt & achieveMax & achieveCov;
                    [~,iMinPlanes] = min(~(planeVec>0));
                    latOptPerf.nPlanesToAchieve(iEcc,iLat) = ExSol.archMat(1,iMinPlanes);
                end
            end
            
            if nSats(iEcc,iSats) <= 80
                % Walker
                load([walkerFolder '\WalkerRgtExSol_Lat_' num2str(latList(iLat))...
                    '_T_' num2str(nSats(iEcc,iSats)) '.mat']);
                walkerPerf.maxPdop(1,iSats,iLat)  = ExSol.maxPdop(ExSol.optF+1,ExSol.optP);
                walkerPerf.intPdop(1,iSats,iLat)  = ExSol.intPdop(ExSol.optF+1,ExSol.optP);
                walkerPerf.coverage(1,iSats,iLat) = ExSol.coverage(ExSol.optF+1,ExSol.optP);
                if isnan(walkerPerf.nSatsToAchieve(iLat))
                    achieveInt = ExSol.intPdop < intTarget;
                    achieveMax = ExSol.maxPdop < maxTarget;
                    achieveCov = ExSol.coverage > covTarget;
                    if any(any(achieveInt & achieveMax & achieveCov))
                        walkerPerf.nSatsToAchieve(iLat) = nSats(1,iSats);
                        planeVec = sum(achieveInt & achieveMax & achieveCov,1);
                        [~,walkerPerf.nPlanesToAchieve(iLat)] = min(~(planeVec>0));
                    end
                end
            end
        end
    end
    % Plot Results for latitude
    figure()% Int PDOP
    semilogy(nSats(1,:),walkerPerf.intPdop(:,:,iLat),'*',...
        nSats.',lat10Perf.intPdop(:,:,iLat).','o',...
        nSats.',latIncPerf.intPdop(:,:,iLat).','x',...
        nSats.',latOptPerf.intPdop(:,:,iLat).','s')
    title(['Integral of PDOP for \phi_0 = ' num2str(latList(iLat))])
    legend('Walker',['Lattice'],'Lattice Inc',['Lattice Optimal'])
    xlabel('# Sats')
    ylabel('$\frac{1}{T}\int^{T}_{0}{PDOPdt}$','interpreter','latex','fontsize',12)
    grid minor
    
    figure()% max PDOP
    semilogy(nSats(1,:),walkerPerf.maxPdop(:,:,iLat),'*',...
        nSats.',lat10Perf.maxPdop(:,:,iLat).','o',...
        nSats.',latIncPerf.maxPdop(:,:,iLat).','x',...
        nSats.',latOptPerf.maxPdop(:,:,iLat).','s')
    title(['Maximum of PDOP for \phi_0 = ' num2str(latList(iLat))])
    legend('Walker',['Lattice'],'Lattice Inc',['Lattice Optimal'])
    xlabel('# Sats')
    ylabel('$\max{PDOP}$','interpreter','latex')
    grid minor
    
    figure()% Coverage
    semilogy(nSats(1,:),walkerPerf.coverage(:,:,iLat),'*',...
        nSats.',lat10Perf.coverage(:,:,iLat).','o',...
        nSats.',latIncPerf.coverage(:,:,iLat).','x',...
        nSats.',latOptPerf.coverage(:,:,iLat).','s')
    title(['Coverage for \phi_0 = ' num2str(latList(iLat))])
    legend('Walker',['Lattice'],'Lattice Inc',['Lattice Optimal'])
    xlabel('# Sats')
    ylabel('Coverage %')
    ylim([0,100])
    grid minor
end

figure()
plot(latList,walkerPerf.nSatsToAchieve(1,:),'*'...
    ,latList,lat10Perf.nSatsToAchieve(1:end,:),'o'...
    ,latList,latIncPerf.nSatsToAchieve(1:end,:),'x'...
    ,latList,latOptPerf.nSatsToAchieve(1:end,:),'s')
title(['Min Sats for: Int PDOP < ' num2str(intTarget) ' & Max PDOP < ' ...
    num2str(maxTarget) ' & Coverage > ' num2str(covTarget)])
xlabel('\phi_0 [°]')
ylabel('# Satellites')
legend('Walker',['Lattice'],'Lattice Inc',['Lattice Optimal'])
grid minor

figure()
plot(latList,walkerPerf.nPlanesToAchieve(1,:),'*'...
    ,latList,lat10Perf.nPlanesToAchieve(1:end,:),'o'...
    ,latList,latIncPerf.nPlanesToAchieve(1:end,:),'o'...
    ,latList,latOptPerf.nPlanesToAchieve(1:end,:),'s')
title(['Min Planes for: Int PDOP < ' num2str(intTarget) ' & Max PDOP < ' ...
    num2str(maxTarget) ' & Coverage > ' num2str(covTarget)])
legend('Walker',['Lattice'],'Lattice Inc',['Lattice Optimal'])
xlabel('\phi_0 [°]')
ylabel('# Planes')
grid minor

%% Find min Planes to Achieve Goal for each nSats
nPlanesToAchieve = nan(4,nCons);
latGs = 20;
for iEcc = 1:length(eccList)
    for iSats = 1:nCons
        % Lattice
        load([lat10Folder '\LatticeExSol_Lat_' num2str(latGs)...
            '_nSats_' num2str(nSats(iEcc,iSats)) '_ecc_' num2str(eccList(iEcc)) '.mat']);
        % Check Goal
        achieveInt = ExSol.intPdop < intTarget;
        achieveMax = ExSol.maxPdop < maxTarget;
        achieveCov = ExSol.coverage > covTarget;
        if any(achieveInt & achieveMax & achieveCov)
            planeVec = achieveInt & achieveMax & achieveCov;
            [~,iMinPlanes] = min(~(planeVec>0));
            nPlanesToAchieve(2,iSats) = ExSol.archMat(1,iMinPlanes);
        end
        
        % Lattice Optimal
        load([latIncFolder '\LatticeExSol_Lat_' num2str(latGs)...
            '_nSats_' num2str(nSats(iEcc,iSats)) '_ecc_' num2str(eccList(iEcc)) '.mat']);
        % Check Goal
        achieveInt = ExSol.intPdop < intTarget;
        achieveMax = ExSol.maxPdop < maxTarget;
        achieveCov = ExSol.coverage > covTarget;
        if any(achieveInt & achieveMax & achieveCov)
            planeVec = achieveInt & achieveMax & achieveCov;
            [~,iMinPlanes] = min(~(planeVec>0));
            nPlanesToAchieve(3,iSats) = ExSol.archMat(1,iMinPlanes);
        end
        
        % Lattice Optimal Inc & Start
        load([latOptFolder '\LatticeExSol_Lat_' num2str(latGs)...
            '_nSats_' num2str(nSats(iEcc,iSats)) '_ecc_' num2str(eccList(iEcc)) '.mat']);
        % Check Goal
        achieveInt = ExSol.intPdop < intTarget;
        achieveMax = ExSol.maxPdop < maxTarget;
        achieveCov = ExSol.coverage > covTarget;
        if any(achieveInt & achieveMax & achieveCov)
            planeVec = achieveInt & achieveMax & achieveCov;
            [~,iMinPlanes] = min(~(planeVec>0));
            nPlanesToAchieve(4,iSats) = ExSol.archMat(1,iMinPlanes);
        end
        
        if nSats(iEcc,iSats) <= 80
            % Walker
            load([walkerFolder '\WalkerRgtExSol_Lat_' num2str(latGs)...
                '_T_' num2str(nSats(iEcc,iSats)) '.mat']);
            achieveInt = ExSol.intPdop < intTarget;
            achieveMax = ExSol.maxPdop < maxTarget;
            achieveCov = ExSol.coverage > covTarget;
            if any(any(achieveInt & achieveMax & achieveCov))
                planeVec = sum(achieveInt & achieveMax & achieveCov,1);
                [~,nPlanesToAchieve(1,iSats)] = min(~(planeVec>0));
            end
        end
    end
end

figure()
plot(nSats(1,:),nPlanesToAchieve(1,:),'*'...
    ,nSats(1,:),nPlanesToAchieve(2,:),'o'...
    ,nSats(1,:),nPlanesToAchieve(3,:),'x'...
    ,nSats(1,:),nPlanesToAchieve(4,:),'s')
title(['Min Planes for: Int PDOP < ' num2str(intTarget) ' & Max PDOP < ' ...
    num2str(maxTarget) ' & Coverage > ' num2str(covTarget)...
    ' \phi_0 = ' num2str(latGs)])
xlabel('# Satellites')
ylabel('# Planes')
legend('Walker',['Lattice'],'Lattice Inc',['Lattice Optimal'],'location','best')
grid minor
