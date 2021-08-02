%% NOTE
% This script has been greatly altered and should be used ONLY to generate
% LFC pareto fronts split by apogee altitude

%% Define Scenarios for Analysis
clear
dataFolder = 'C:\Users\User\Google Drive\Master''s Degree\Lattice Optimization Data\Final Results';
folderList = {...
    'C:\Users\User\Google Drive\Master''s Degree\Walker Optimization Data\Previous Optimization Runs\Latticified Walker multilat multiInc collision';...
    'C:\Users\User\Google Drive\Master''s Degree\Lattice Optimization Data\Previous Runs\LatticeDef v1';...
    'C:\Users\User\Google Drive\Master''s Degree\Lattice Optimization Data\GA Standard\Previous Runs\Version 6 - definitive';...
    };
markerList = {...
    'o';...
    'v';...
    'x';...
    };
nameList = {...
    'Circular, h_a = h_c';...
    'Elliptical, h_a = 900 km';...
    'Elliptical, h_a = 1000 km';...
    };
nScenarios = numel(folderList);

%% Set Performance Goals
intTarget = 5;
maxTarget = 1e10;
covTarget = 95;
p90Target = 1e10;

%% Analyze Scenarios
close all
nLeg = 0;
for iScenario = 2 %1:nScenarios
    load([folderList{iScenario} '\OptParams.mat']);
    if exist('OptParams','var')
        % adjust for different data structure
        hAList = 1; % No separation of results by hA
        PropParams.maxPdop = OptParams.maxPdop;
        PropParams.timeVec = OptParams.timeVec;
        PropParams.elevMin = OptParams.elevMin;
        PropParams.relTol = OptParams.relTol;
        PropParams.absTol = OptParams.absTol;
        PropParams.datafolder = OptParams.datafolder;
        PropParams.delLat = OptParams.delLat;
        clear OptParams;
        gaFlag = true;
    else
        gaFlag = false;
    end
    nSats = minSats:maxSats;
    maxPdop = nan(numel(hAList),numel(nSats));
    intPdop = nan(numel(hAList),numel(nSats));
    coverage = nan(numel(hAList),numel(nSats));
    p90 = nan(numel(hAList),numel(nSats));
%     latList = 50;
%     hAList = 1000;
    for iLat = 1:numel(latList)
        paretoSats = mat2cell(inf(numel(hAList),1),ones(1,numel(hAList)));
        paretoPlanes = mat2cell(inf(numel(hAList),1),ones(1,numel(hAList)));
        if ~gaFlag
            for iHA = 1:numel(hAList)
                for iSats = 1:numel(nSats)
                    load([folderList{iScenario} '\LatticeExSol_Lat_'...
                        num2str(latList(iLat)) '_nSats_' num2str(nSats(iSats)) ...
                        '_hA_' num2str(hAList(iHA)) '.mat']);
                    
                    achieveInt = ExSol.intPdop < intTarget;
                    achieveMax = ExSol.maxPdop < maxTarget;
                    achieveCov = ExSol.coverage > covTarget;
                    achievep90 = ExSol.p90 < p90Target;
                    if any(achieveInt & achieveMax & achieveCov & achievep90)
                        planeVec = achieveInt & achieveMax & achieveCov & achievep90;
                        [~,iMinPlanes] = min(~(planeVec>0));
                        nPlanesToAchieve = ExSol.archMat(1,iMinPlanes);
                        if nPlanesToAchieve < min(paretoPlanes{iHA})
                            if min(paretoPlanes{iHA}) < inf
                                paretoSats{iHA} = [paretoSats{iHA}, nSats(iSats)];
                                paretoPlanes{iHA} = [paretoPlanes{iHA}, nPlanesToAchieve];
                            else
                                paretoSats{iHA} = [nSats(iSats)];
                                paretoPlanes{iHA} = [nPlanesToAchieve];
                            end
                        end
                    end
                    maxPdop(iHA,iSats) = ExSol.maxPdop(ExSol.iOpt);
                    intPdop(iHA,iSats) = ExSol.intPdop(ExSol.iOpt);
                    coverage(iHA,iSats) = ExSol.coverage(ExSol.iOpt);
                    p90(iHA,iSats) = ExSol.p90(ExSol.iOpt);
                end
            end
        else % Genetic Algorithm
            iHA = 1;
            for iSats = 1:numel(nSats)
                load([folderList{iScenario} '\LatticeGaSol_Lat_'...
                    num2str(latList(iLat)) '_nSats_' num2str(nSats(iSats)) '.mat']);
                
                achieveInt = mean(GaSol.intPdop,1) < intTarget;
                achieveMax = mean(GaSol.maxPdop,1) < maxTarget;
                achieveCov = mean(GaSol.coverage,1) > covTarget;
                achievep90 = mean(GaSol.p90,1) < p90Target;
                if any(achieveInt & achieveMax & achieveCov & achievep90)
                    planeVec = achieveInt & achieveMax & achieveCov & achievep90;
                    [~,iMinPlanes] = min(~(planeVec>0));
                    nPlanesToAchieve = GaSol.archMat(1,iMinPlanes);
                    if nPlanesToAchieve < min(paretoPlanes{iHA})
                        if min(paretoPlanes{iHA}) < inf
                            paretoSats{iHA} = [paretoSats{iHA}, nSats(iSats)];
                            paretoPlanes{iHA} = [paretoPlanes{iHA}, nPlanesToAchieve];
                        else
                            paretoSats{iHA} = [nSats(iSats)];
                            paretoPlanes{iHA} = [nPlanesToAchieve];
                        end
                    end
                end
                maxPdop(iHA,iSats) = mean(GaSol.maxPdop(:,GaSol.iOpt),1);
                intPdop(iHA,iSats) = mean(GaSol.intPdop(:,GaSol.iOpt),1);
                coverage(iHA,iSats) = mean(GaSol.coverage(:,GaSol.iOpt),1);
                p90(iHA,iSats) = mean(GaSol.p90(:,GaSol.iOpt),1);
            end
        end
        
%         % Plot Results for Latitude
%         figure(iLat*10 + 1) % Int PDOP
%         hold on
%         semilogy(nSats,intPdop,markerList{iScenario})
%         title(['Integral of PDOP for \phi_0 = ' num2str(latList(iLat))])
%         xlabel('# Sats')
%         ylabel('$\frac{1}{T}\int^{T}_{0}{PDOPdt}$','interpreter','latex','fontsize',12)
%         grid on
%         set(gca, 'YScale', 'log')
%         leg = get(gca,'Legend');
%         if isempty(leg)
%             leg = legend(nameList{iScenario});
%         else
%             leg.String(nLeg+1:nLeg + numel(hAList)) = nameList(iScenario);
%         end
%         hold off
        
%         figure(iLat*10 + 2) % 90th Percentile PDOP
%         hold on
%         semilogy(nSats,p90,markerList{iScenario})
%         title(['90th Percentile of PDOP for \phi_0 = ' num2str(latList(iLat))])
%         xlabel('# Sats')
%         ylabel('$PDOP$','interpreter','latex','fontsize',12)
%         grid on
%         set(gca, 'YScale', 'log')
%         leg = get(gca,'Legend');
%         if isempty(leg)
%             leg = legend(nameList{iScenario});
%         else
%             leg.String(nLeg+1:nLeg + numel(hAList)) = nameList(iScenario);
%         end
%         hold off
        
%         figure(iLat*10 + 3) % Coverage
%         hold on
%         plot(nSats,coverage,markerList{iScenario})
%         title(['Coverage for \phi_0 = ' num2str(latList(iLat))])
%         xlabel('# Sats')
%         ylabel('Coverage %')
%         ylim([0,100])
%         grid on
%         leg = get(gca,'Legend');
%         if isempty(leg)
%             leg = legend(nameList{iScenario});
%         else
%             leg.String(nLeg+1:nLeg + numel(hAList)) = nameList(iScenario);
%         end
%         hold off
        
        figure(iLat*10 + 4) % Pareto Frontier
        hold on
        for iHA = 1:numel(hAList)
            plot(paretoSats{iHA},paretoPlanes{iHA},markerList{iHA}...
                ,'linewidth',1.5,'markersize',10)
        end
%         title(['Pareto Frontier \phi_0 = ' num2str(latList(iLat))])
        xlabel('$N_s$','interpreter','latex','fontsize',14)
        ylabel('$N_o$','interpreter','latex','fontsize',14)
        xlim([50,80])
        ylim([0,65])
        grid on
        legend(nameList)
%         leg = get(gca,'Legend');
%         if isempty(leg)
%             leg = legend(nameList{iHA});
%         else
%             leg.String(nLeg) = nameList(iHA);
%         end
        hold off
        
        paretoSet = struct();
        paretoSet.latEm = latList(iLat);
        paretoSet.nSats = paretoSats;
        paretoSet.nPlanes = paretoPlanes;
        paretoSet.intTarget = intTarget;
        paretoSet.maxTarget = maxTarget;
        paretoSet.p90Target = p90Target;
        paretoSet.covTarget = covTarget;
        
        save([dataFolder '\ParetoList_Lat_' num2str(paretoSet.latEm) ...
            '_LFC Ex_hASplit.mat'],'paretoSet');
        
    end
    nLeg = nLeg + numel(hAList);
end