%% Define Scenarios for Analysis
folderList = {'C:\Users\User\Dropbox\Walker Optimization Data\Previous Optimization Runs\Latticified Walker',...
    'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\Apogee Height, delta inc',...
    'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\Apogee Height, opt inc(pre-fix)',...
    'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\Apogee Height x3, opt inc'};
markerList = {'*','o','x','s'};
nameList = {'Walker','Lattice Basic','Lattice Optimal 1','Lattice Optimal 2'};
nScenarios = numel(folderList);

%% Set Performance Goals
intTarget = 10;
maxTarget = 1000;
covTarget = 99;
p90Target = 5;

%% Analyze Scenarios
close all
nLeg = 0;
for iScenario = 1:nScenarios
    load([folderList{iScenario} '\OptParams.mat']);
    hAList = 0;
    nSats = minSats:maxSats;
    maxPdop = nan(numel(hAList),numel(nSats));
    intPdop = nan(numel(hAList),numel(nSats));
    coverage = nan(numel(hAList),numel(nSats));
    p90 = nan(numel(hAList),numel(nSats));
    for iLat = 1:numel(latList)
        paretoSats = mat2cell(inf(numel(hAList),1),ones(1,numel(hAList)));
        paretoPlanes = mat2cell(inf(numel(hAList),1),ones(1,numel(hAList)));
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
                        paretoSats{iHA} = [paretoSats{iHA}, nSats(iSats)];
                        paretoPlanes{iHA} = [paretoPlanes{iHA}, nPlanesToAchieve];
                    end
                end
                maxPdop(iHA,iSats) = ExSol.maxPdop(ExSol.iOpt);
                intPdop(iHA,iSats) = ExSol.intPdop(ExSol.iOpt);
                coverage(iHA,iSats) = ExSol.coverage(ExSol.iOpt);
                p90(iHA,iSats) = ExSol.p90(ExSol.iOpt);
            end
        end
        % Plot Results for Latitude
        figure(iLat*10 + 1) % Int PDOP
        hold on
        semilogy(nSats,intPdop,markerList{iScenario})
        title(['Integral of PDOP for \phi_0 = ' num2str(latList(iLat))])
        xlabel('# Sats')
        ylabel('$\frac{1}{T}\int^{T}_{0}{PDOPdt}$','interpreter','latex','fontsize',12)
        grid on
        set(gca, 'YScale', 'log')
        leg = get(gca,'Legend');
        if isempty(leg)
            leg = legend(nameList{iScenario});
        else
            leg.String(nLeg+1:nLeg + numel(hAList)) = nameList(iScenario);
        end
        hold off
        
        figure(iLat*10 + 2) % 90th Percentile PDOP
        hold on
        semilogy(nSats,p90,markerList{iScenario})
        title(['90th Percentile of PDOP for \phi_0 = ' num2str(latList(iLat))])
        xlabel('# Sats')
        ylabel('$PDOP$','interpreter','latex','fontsize',12)
        grid on
        set(gca, 'YScale', 'log')
        leg = get(gca,'Legend');
        if isempty(leg)
            leg = legend(nameList{iScenario});
        else
            leg.String(nLeg+1:nLeg + numel(hAList)) = nameList(iScenario);
        end
        hold off
        
        figure(iLat*10 + 3) % Coverage
        hold on
        plot(nSats,coverage,markerList{iScenario})
        title(['Coverage for \phi_0 = ' num2str(latList(iLat))])
        xlabel('# Sats')
        ylabel('Coverage %')
        ylim([0,100])
        grid on
        leg = get(gca,'Legend');
        if isempty(leg)
            leg = legend(nameList{iScenario});
        else
            leg.String(nLeg+1:nLeg + numel(hAList)) = nameList(iScenario);
        end
        hold off
        
        figure(iLat*10 + 4) % Pareto Frontier
        hold on
        plot(paretoSats{:},paretoPlanes{:},['--' markerList{iScenario}]...
        ,'linewidth',1.5,'markersize',10)
        title(['Pareto Frontier \phi_0 = ' num2str(latList(iLat))])
        xlabel('# Sats')
        ylabel('# Planes')
        grid on
        leg = get(gca,'Legend');
        if isempty(leg)
            leg = legend(nameList{iScenario});
        else
            leg.String(nLeg+1:nLeg + numel(hAList)) = nameList(iScenario);
        end
        hold off
    end
    nLeg = nLeg + numel(hAList);
end