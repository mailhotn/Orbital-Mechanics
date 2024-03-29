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
    'Walker exhaustive search';...
    'LFC exhaustive search';...
    'LFC genetic algorithm';...
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
for iScenario = 1:nScenarios
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
        paretoSats = inf;
        paretoPlanes = inf;
        paretoHAs = [];
        if ~gaFlag
            for iSats = 1:numel(nSats)
                for iHA = 1:numel(hAList)
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
                        if nPlanesToAchieve < min(paretoPlanes)
                            if min(paretoPlanes) < inf
                                if paretoSats(end) < nSats(iSats)
                                    paretoSats = [paretoSats, nSats(iSats)];
                                    paretoPlanes = [paretoPlanes, nPlanesToAchieve];
                                    paretoHAs = [paretoHAs, hAList(iHA)];
                                else
                                    paretoSats(end) = nSats(iSats);
                                    paretoPlanes(end) = nPlanesToAchieve;
                                    paretoHAs(end) = hAList(iHA);
                                end
                            else
                                paretoSats = [nSats(iSats)];
                                paretoPlanes = [nPlanesToAchieve];
                                paretoHAs = hAList(iHA);
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
                    if nPlanesToAchieve < min(paretoPlanes)
                        if min(paretoPlanes) < inf
                            paretoSats = [paretoSats, nSats(iSats)];
                            paretoPlanes = [paretoPlanes, nPlanesToAchieve];
                        else
                            paretoSats = [nSats(iSats)];
                            paretoPlanes = [nPlanesToAchieve];
                        end
                    end
                end
                maxPdop(iHA,iSats) = mean(GaSol.maxPdop(:,GaSol.iOpt),1);
                intPdop(iHA,iSats) = mean(GaSol.intPdop(:,GaSol.iOpt),1);
                coverage(iHA,iSats) = mean(GaSol.coverage(:,GaSol.iOpt),1);
                p90(iHA,iSats) = mean(GaSol.p90(:,GaSol.iOpt),1);
            end
        end
        
        % Plot Results for Latitude
        figure(iLat*10 + 1) % Int PDOP
        hold on
        semilogy(nSats,min(intPdop,[],1),markerList{iScenario})
        if iScenario == 3
            semilogy([minSats,maxSats],[intTarget,intTarget],'b','linewidth',1.5)
        end
%         title(['Integral of PDOP for \phi_0 = ' num2str(latList(iLat))])
        xlabel('$N_s$','interpreter','latex','fontsize',14)
%         ylabel('$\frac{1}{T}\int^{T}_{0}{PDOPdt}$','interpreter','latex','fontsize',12)
        ylabel('$\mathcal{J}$','interpreter','latex','fontsize',14)
        grid on
        set(gca, 'YScale', 'log')
        legend(nameList)
%         leg = get(gca,'Legend');
%         if isempty(leg)
%             leg = legend(nameList{iScenario});
%         else
%             leg.String(nLeg) = nameList(iScenario);
%         end
        hold off
        
%         figure(iLat*10 + 2) % 90th Percentile PDOP
%         hold on
%         semilogy(nSats,min(p90,[],1),markerList{iScenario})
%         title(['90th Percentile of PDOP for \phi_0 = ' num2str(latList(iLat))])
%         xlabel('# Sats')
%         ylabel('$PDOP$','interpreter','latex','fontsize',12)
%         grid on
%         set(gca, 'YScale', 'log')
%         leg = get(gca,'Legend');
%         if isempty(leg)
%             leg = legend(nameList{iScenario});
%         else
%             leg.String(nLeg+1) = nameList(iScenario);
%         end
%         hold off
%         
%         figure(iLat*10 + 3) % Coverage
%         hold on
%         plot(nSats,min(coverage,[],1),markerList{iScenario})
%         title(['Coverage for \phi_0 = ' num2str(latList(iLat))])
%         xlabel('# Sats')
%         ylabel('Coverage %')
%         ylim([0,100])
%         grid on
%         leg = get(gca,'Legend');
%         if isempty(leg)
%             leg = legend(nameList{iScenario});
%         else
%             leg.String(nLeg+1) = nameList(iScenario);
%         end
%         hold off
        
        figure(iLat*10 + 4) % Pareto Frontier
        hold on
            plot(paretoSats,paretoPlanes,markerList{iScenario}...
                ,'linewidth',1.5,'markersize',10)
%         title(['Pareto Frontier \phi_0 = ' num2str(latList(iLat))])
        xlabel('$N_s$','interpreter','latex','fontsize',14)
        ylabel('$N_o$','interpreter','latex','fontsize',14)
        xlim([45,80])
        ylim([0,65])
        grid on
        legend(nameList)
%         leg = get(gca,'Legend');
%         if isempty(leg)
%             leg = legend(nameList{iScenario});
%         else
%             leg.String(nLeg+1) = nameList(iScenario);
%         end
        hold off
        
        paretoSet = struct();
        paretoSet.latEm = latList(iLat);
        paretoSet.nSats = paretoSats;
        paretoSet.nPlanes = paretoPlanes;
        paretoSet.hAs = paretoHAs;
        paretoSet.intTarget = intTarget;
        paretoSet.maxTarget = maxTarget;
        paretoSet.p90Target = p90Target;
        paretoSet.covTarget = covTarget;
        
        save([dataFolder '\ParetoList_Lat_' num2str(paretoSet.latEm) '_'...
            nameList{iScenario} '.mat'],'paretoSet');
        
    end
    nLeg = nLeg + numel(hAList);
end