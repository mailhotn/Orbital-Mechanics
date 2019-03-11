%% Load Optimization Parameters

datafolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\Lattice Version 2';
walkerFolder = ['C:\Users\User\Dropbox\Walker Optimization Data'...
    '\Previous Optimization Runs\Walker RGT Ex Search delta inc 10'];

load ([datafolder '\OptParams.mat']);
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

latticeMaxPdop  = nan(nEcc,nCons,nLats);
latticeIntPdop  = nan(nEcc,nCons,nLats);
latticeCoverage = nan(nEcc,nCons,nLats);
walkerMaxPdop   = nan(1,nCons,nLats);
walkerIntPdop   = nan(1,nCons,nLats);
walkerCoverage  = nan(1,nCons,nLats);
% inc       = nan(nLats,nCons);
% alt       = nan(nLats,nCons);
% latticeNSats    = repmat(minSats:maxSats,nLats,1);
nSats           = repmat(minSats:maxSats,nEcc,1);

%% Performance Goals
intTarget = 5;
maxTarget = 20;
covTarget = 99.5;
nSatsToAchieve = nan(nEcc+1,nLats);
nPlanesToAchieve = nan(nEcc+1,nLats);

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
            % Lattice
            load([datafolder '\LatticeExSol_Lat_' num2str(latList(iLat))...
                '_nSats_' num2str(nSats(iEcc,iSats)) '_ecc_' num2str(eccList(iEcc)) '.mat']);
            latticeMaxPdop(iEcc,iSats,iLat) = ExSol.maxPdop(ExSol.iOpt);
            latticeIntPdop(iEcc,iSats,iLat) = ExSol.intPdop(ExSol.iOpt);
            latticeCoverage(iEcc,iSats,iLat) = ExSol.coverage(ExSol.iOpt);
            % Check Goal
            if isnan(nSatsToAchieve(iEcc+1,iLat))
                achieveInt = ExSol.intPdop < intTarget;
                achieveMax = ExSol.maxPdop < maxTarget;
                achieveCov = ExSol.coverage > covTarget;
                if any(achieveInt & achieveMax & achieveCov)
                    nSatsToAchieve(iEcc+1, iLat) = nSats(iEcc,iSats);
                    planeVec = achieveInt & achieveMax & achieveCov;
                    [~,iMinPlanes] = min(~(planeVec>0));
                    nPlanesToAchieve(iEcc+1,iLat) = ExSol.archMat(1,iMinPlanes);
                end
            end
            if nSats(iEcc,iSats) <= 80
                % Walker
                load([walkerFolder '\WalkerRgtExSol_Lat_' num2str(latList(iLat))...
                    '_T_' num2str(nSats(iEcc,iSats)) '.mat']);
                walkerMaxPdop(1,iSats,iLat)  = ExSol.maxPdop(ExSol.optF+1,ExSol.optP);
                walkerIntPdop(1,iSats,iLat)  = ExSol.intPdop(ExSol.optF+1,ExSol.optP);
                walkerCoverage(1,iSats,iLat) = ExSol.coverage(ExSol.optF+1,ExSol.optP);
                if isnan(nSatsToAchieve(1,iLat))
                    achieveInt = ExSol.intPdop < intTarget;
                    achieveMax = ExSol.maxPdop < maxTarget;
                    achieveCov = ExSol.coverage > covTarget;
                    if any(any(achieveInt & achieveMax & achieveCov))
                        nSatsToAchieve(1, iLat) = nSats(1,iSats);
                        planeVec = sum(achieveInt & achieveMax & achieveCov,1);
                        [~,nPlanesToAchieve(1,iLat)] = min(~(planeVec>0));
                    end
                end
            end
        end
    end
    % Plot Results for latitude
    figure()% Int PDOP
    semilogy(nSats(1,:),walkerIntPdop(:,:,iLat),'*',nSats.',latticeIntPdop(:,:,iLat).','o')
    title(['Integral of PDOP for \phi_0 = ' num2str(latList(iLat))])
    legend('Walker',['Lattice e = ' num2str(eccList(1))],['Lattice e = ' num2str(eccList(2))],...
        ['Lattice e = ' num2str(eccList(3))])
    xlabel('# Sats')
    ylabel('$\frac{1}{T}\int^{T}_{0}{PDOPdt}$','interpreter','latex','fontsize',12)
    grid minor
    
    figure()% Max PDOP
    semilogy(nSats(1,:),walkerMaxPdop(:,:,iLat),'*',nSats.',latticeMaxPdop(:,:,iLat).','o')
    title(['Maximum of PDOP for \phi_0 = ' num2str(latList(iLat))])
    legend('Walker',['Lattice e = ' num2str(eccList(1))],['Lattice e = ' num2str(eccList(2))],...
        ['Lattice e = ' num2str(eccList(3))])
    xlabel('# Sats')
    ylabel('$\max{PDOP}$','interpreter','latex')
    grid minor
    
    figure()% Coverage
    plot(nSats(1,:),walkerCoverage(:,:,iLat),'*',nSats.',latticeCoverage(:,:,iLat).','o')
    title(['Coverage for \phi_0 = ' num2str(latList(iLat))])
    legend('Walker',['Lattice e = ' num2str(eccList(1))],['Lattice e = ' num2str(eccList(2))],...
        ['Lattice e = ' num2str(eccList(3))])
    xlabel('# Sats')
    ylabel('Coverage %')
    ylim([0,100])
    grid minor
end

figure()
plot(latList,nSatsToAchieve(1,:),'*'...
    ,latList,nSatsToAchieve(2:end,:),'o')
title(['Min Sats for: Int PDOP < ' num2str(intTarget) ' & Max PDOP < ' ...
    num2str(maxTarget) ' & Coverage > ' num2str(covTarget)])
xlabel('\phi_0 [°]')
ylabel('# Satellites')
legend('Walker',['Lattice e = ' num2str(eccList(1))],['Lattice e = ' num2str(eccList(2))],...
    ['Lattice e = ' num2str(eccList(3))],['Lattice e = ' num2str(eccList(4))],'location','best')
grid minor

figure()
plot(latList,nPlanesToAchieve(1,:),'*'...
    ,latList,nPlanesToAchieve(2:end,:),'o')
title(['Min Planes for: Int PDOP < ' num2str(intTarget) ' & Max PDOP < ' ...
    num2str(maxTarget) ' & Coverage > ' num2str(covTarget)])
xlabel('\phi_0 [°]')
ylabel('# Planes')
legend('Walker',['Lattice e = ' num2str(eccList(1))],['Lattice e = ' num2str(eccList(2))],...
    ['Lattice e = ' num2str(eccList(3))],['Lattice e = ' num2str(eccList(4))],'location','best')
grid minor

%% Find min Planes to Achieve Goal for each nSats
nPlanesToAchieve = nan(nEcc+1,nCons);
latGs = 50;
for iEcc = 1:length(eccList)
    for iSats = 1:nCons
        % Lattice
        load([datafolder '\LatticeExSol_Lat_' num2str(latGs)...
            '_nSats_' num2str(nSats(iEcc,iSats)) '_ecc_' num2str(eccList(iEcc)) '.mat']);
        % Check Goal
        achieveInt = ExSol.intPdop < intTarget;
        achieveMax = ExSol.maxPdop < maxTarget;
        achieveCov = ExSol.coverage > covTarget;
        if any(achieveInt & achieveMax & achieveCov)
            planeVec = achieveInt & achieveMax & achieveCov;
            [~,iMinPlanes] = min(~(planeVec>0));
            nPlanesToAchieve(iEcc+1,iSats) = ExSol.archMat(1,iMinPlanes);
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
    ,nSats(1,:),nPlanesToAchieve(2:end,:),'o')
title(['Min Planes for: Int PDOP < ' num2str(intTarget) ' & Max PDOP < ' ...
    num2str(maxTarget) ' & Coverage > ' num2str(covTarget)...
    ' \phi_0 = ' num2str(latGs)])
xlabel('# Satellites')
ylabel('# Planes')
legend('Walker',['Lattice e = ' num2str(eccList(1))],['Lattice e = ' num2str(eccList(2))],...
    ['Lattice e = ' num2str(eccList(3))],['Lattice e = ' num2str(eccList(4))],'location','best')
grid minor
