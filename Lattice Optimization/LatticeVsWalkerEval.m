%% Load Optimization Parameters


datafolder = ['C:\Users\User\Dropbox\Graduate Research Day 2019'...
    '\Optimization Data\Lattice Alt lat 30'];
walkerFolder = ['C:\Users\User\Dropbox\Graduate Research Day 2019'...
    '\Optimization Data\Walker lat 30'];


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
nCons = maxSats - minSats + 1;
nHA = length(hAList);

latticeMaxPdop  = nan(nHA,nCons);
latticeIntPdop  = nan(nHA,nCons);
latticeCoverage = nan(nHA,nCons);
walkerMaxPdop   = nan(1,nCons);
walkerIntPdop   = nan(1,nCons);
walkerCoverage  = nan(1,nCons);
% inc       = nan(nLats,nCons);
% alt       = nan(nLats,nCons);
% latticeNSats    = repmat(minSats:maxSats,nLats,1);
nSats           = repmat(minSats:maxSats,nHA,1);

%% Performance Goals
intTarget = 3;
maxTarget = 10000000;
covTarget = 10;
nSatsToAchieve = nan(nHA+1,1);
nPlanesToAchieve = nan(nHA+1,1);

%% Check Performance
%     latticeMaxPdop  = nan(nEcc,nCons);
%     latticeIntPdop  = nan(nEcc,nCons);
%     latticeCoverage = nan(nEcc,nCons);
%     walkerMaxPdop   = nan(1,nCons);
%     walkerIntPdop   = nan(1,nCons);
%     walkerCoverage  = nan(1,nCons);
for iHA = 1:length(hAList)
    for iSats = 1:nCons
        % Lattice
        load([datafolder '\LatticeExSol_Lat_' num2str(latEm)...
            '_nSats_' num2str(nSats(iHA,iSats)) '_hA_' num2str(hAList(iHA)) '.mat']);
        latticeMaxPdop(iHA,iSats) = ExSol.maxPdop(ExSol.iOpt);
        latticeIntPdop(iHA,iSats) = ExSol.intPdop(ExSol.iOpt);
        latticeCoverage(iHA,iSats) = ExSol.coverage(ExSol.iOpt);
        % Check Goal
        if isnan(nSatsToAchieve(iHA+1))
            achieveInt = ExSol.intPdop < intTarget;
            achieveMax = ExSol.maxPdop < maxTarget;
            achieveCov = ExSol.coverage > covTarget;
            if any(achieveInt & achieveMax & achieveCov)
                nSatsToAchieve(iHA+1) = nSats(iHA,iSats);
                planeVec = achieveInt & achieveMax & achieveCov;
                [~,iMinPlanes] = min(~(planeVec>0));
                nPlanesToAchieve(iHA+1) = ExSol.archMat(1,iMinPlanes);
            end
        end
        if nSats(iHA,iSats) <= 80
            % Walker
            load([walkerFolder '\WalkerRgtExSol_Lat_' num2str(latEm)...
                '_T_' num2str(nSats(iHA,iSats)) '.mat']);
            walkerMaxPdop(1,iSats)  = ExSol.maxPdop(ExSol.optF+1,ExSol.optP);
            walkerIntPdop(1,iSats)  = ExSol.intPdop(ExSol.optF+1,ExSol.optP);
            walkerCoverage(1,iSats) = ExSol.coverage(ExSol.optF+1,ExSol.optP);
            if isnan(nSatsToAchieve(1))
                achieveInt = ExSol.intPdop < intTarget;
                achieveMax = ExSol.maxPdop < maxTarget;
                achieveCov = ExSol.coverage > covTarget;
                if any(any(achieveInt & achieveMax & achieveCov))
                    nSatsToAchieve(1) = nSats(1,iSats);
                    planeVec = sum(achieveInt & achieveMax & achieveCov,1);
                    [~,nPlanesToAchieve(1)] = min(~(planeVec>0));
                end
            end
        end
    end
end
% Plot Results for latitude
figure()% Int PDOP
semilogy(nSats(1,:),walkerIntPdop,'*',nSats.',latticeIntPdop.','o')
% title(['Integral of PDOP for \phi_0 = ' num2str(latEm)])
legend('Walker','Lattice Circular',['Lattice h_a = ' num2str(hAList(2))],...
    ['Lattice h_a = ' num2str(hAList(3))])
xlabel('# Satellites')
ylabel('$\frac{1}{T}\int^{T}_{0}{PDOP(t)dt}$','interpreter','latex','fontsize',12)
grid on
% 
% figure()% Max PDOP
% semilogy(nSats(1,:),walkerMaxPdop,'*',nSats.',latticeMaxPdop.','o')
% title(['Maximum of PDOP for \phi_0 = ' num2str(latEm)])
% legend('Walker','Lattice Circular',['Lattice h_a = ' num2str(hAList(2))],...
%     ['Lattice h_a = ' num2str(hAList(3))])
% xlabel('# Satellites')
% ylabel('$\max{PDOP}$','interpreter','latex')
% grid on

% figure()% Coverage
% plot(nSats(1,:),walkerCoverage,'*',nSats.',latticeCoverage.','o')
% title(['Coverage for \phi_0 = ' num2str(latEm)])
% legend('Walker','Lattice Circular',['Lattice h_a = ' num2str(hAList(2))],...
%     ['Lattice h_a = ' num2str(hAList(3))])
% xlabel('# Satellites')
% ylabel('Coverage %')
% ylim([0,100])
% grid on

%% Find min Planes to Achieve Goal for each nSats

nPlanesToAchieve = nan(nHA+1,nCons);
for iHA = 1:length(hAList)

    for iSats = 1:nCons
        % Lattice
        load([datafolder '\LatticeExSol_Lat_' num2str(latEm)...
            '_nSats_' num2str(nSats(iHA,iSats)) '_hA_' num2str(hAList(iHA)) '.mat']);
        % Check Goal
        achieveInt = ExSol.intPdop < intTarget;
        achieveMax = ExSol.maxPdop < maxTarget;
        achieveCov = ExSol.coverage > covTarget;
        if any(achieveInt & achieveMax & achieveCov)
            planeVec = achieveInt & achieveMax & achieveCov;
            [~,iMinPlanes] = min(~(planeVec>0));
            nPlanesToAchieve(iHA+1,iSats) = ExSol.archMat(1,iMinPlanes);
        end
        
        if nSats(iHA,iSats) <= 80
            % Walker
            load([walkerFolder '\WalkerRgtExSol_Lat_' num2str(latEm)...
                '_T_' num2str(nSats(iHA,iSats)) '.mat']);
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
for iHA = 1:length(hAList)+1
    for iSats = 2:nCons % Pareto
        if nPlanesToAchieve(iHA,iSats) >= min(nPlanesToAchieve(iHA,1:iSats-1))
            nPlanesToAchieve(iHA,iSats) = nan;
        end
    end
end
figure()
plot(nSats(1,~isnan(nPlanesToAchieve(1,:)))...
    ,nPlanesToAchieve(1,~isnan(nPlanesToAchieve(1,:))),'-s'...
    ,nSats(1,~isnan(nPlanesToAchieve(2,:)))...
    ,nPlanesToAchieve(2,~isnan(nPlanesToAchieve(2,:))),'--x'...
    ,nSats(1,~isnan(nPlanesToAchieve(3,:)))...
    ,nPlanesToAchieve(3,~isnan(nPlanesToAchieve(3,:))),'--x'...
    ,nSats(1,~isnan(nPlanesToAchieve(4,:)))...
    ,nPlanesToAchieve(4,~isnan(nPlanesToAchieve(4,:))),'--x',...
    'linewidth',1.5,'markersize',10)
% title(['Min Planes for: Int PDOP < ' num2str(intTarget) ' & Max PDOP < ' ...
%     num2str(maxTarget) ' & Coverage > ' num2str(covTarget)...
%     ' \phi_0 = ' num2str(latEm)])
xlabel('# Satellites')
ylabel('# Planes')
legend('Walker','Lattice Circular',['Lattice h_a = ' num2str(hAList(2))],...
    ['Lattice h_a = ' num2str(hAList(3))])
grid on
