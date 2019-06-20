clear
close all
optDataFolder = 'C:\Users\User\Dropbox\Walker Optimization Data\Previous Optimization Runs\Walker lat30-60 dInc test';
nomDataFolder = 'C:\Users\User\Dropbox\Walker Optimization Data\Previous Optimization Runs\Latticified Walker';
%% Set Performance Goals
intTarget = 10;
maxTarget = 10;
covTarget = 99;
p90Target = 3;
%% Optimal Walker
load([optDataFolder '\OptParams.mat']);

nSats = minSats:maxSats;
maxPdop = nan(numel(hAList),numel(nSats));
intPdop = nan(numel(hAList),numel(nSats));
coverage = nan(numel(hAList),numel(nSats));
p90 = nan(numel(hAList),numel(nSats));
for iLat = 1:numel(latList)
    paretoSats = inf;
    paretoPlanes = inf;
    paretoDelInc = inf;
    paretoPerf = inf;
    for iSats = 1:numel(nSats)
        for iHA = 1:numel(hAList)
            load([optDataFolder '\LatticeExSol_Lat_'...
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
                if (nPlanesToAchieve < min(paretoPlanes)) && ...
                        (nSats(iSats) > max(paretoSats))
                    
                    paretoSats = [paretoSats, nSats(iSats)];
                    paretoPlanes = [paretoPlanes, nPlanesToAchieve];
                    paretoDelInc = [paretoDelInc, hAList(iHA)];
                    paretoPerf = [paretoPerf, ExSol.p90(iMinPlanes)];
                elseif (nPlanesToAchieve < min(paretoPlanes)) && ...
                        (nSats(iSats) <= max(paretoSats))
                    
                    paretoSats(end) = nSats(iSats);
                    paretoPlanes(end) = nPlanesToAchieve;
                    paretoDelInc(end) = hAList(iHA);
                    paretoPerf(end) = ExSol.p90(iMinPlanes);
                elseif (nPlanesToAchieve == min(paretoPlanes)) &&...
                        (ExSol.p90(iMinPlanes) < paretoPerf(end)) &&...
                        (nSats(iSats) <= paretoSats(end))
                    
                    paretoSats(end) = nSats(iSats);
                    paretoPlanes(end) = nPlanesToAchieve;
                    paretoDelInc(end) = hAList(iHA);
                    paretoPerf(end) = ExSol.p90(iMinPlanes);
                end
            end
            maxPdop(iHA,iSats) = ExSol.maxPdop(ExSol.iOpt);
            intPdop(iHA,iSats) = ExSol.intPdop(ExSol.iOpt);
            coverage(iHA,iSats) = ExSol.coverage(ExSol.iOpt);
            p90(iHA,iSats) = ExSol.p90(ExSol.iOpt);
        end
    end
    figure(iLat*10+1)
    plot(paretoSats,paretoPlanes,'--o')
    title(['Pareto Front \phi_0 = ' num2str(latList(iLat))])
    text(paretoSats+0.5,paretoPlanes+0.5,num2str(paretoDelInc(:)))
    xlabel('# Sats')
    ylabel('# Planes')
    xlim([50,80])
    ylim([0,60])
    grid on
end
%% Nominal Walker
load([nomDataFolder '\OptParams.mat']);
nSats = minSats:maxSats;
maxPdop = nan(numel(hAList),numel(nSats));
intPdop = nan(numel(hAList),numel(nSats));
coverage = nan(numel(hAList),numel(nSats));
p90 = nan(numel(hAList),numel(nSats));
latList = 30:10:60;
for iLat = 1:numel(latList)
    paretoSats = mat2cell(inf(numel(hAList),1),ones(1,numel(hAList)));
    paretoPlanes = mat2cell(inf(numel(hAList),1),ones(1,numel(hAList)));
    for iHA = 1:numel(hAList)
        for iSats = 1:numel(nSats)
            load([nomDataFolder '\LatticeExSol_Lat_'...
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
    figure(iLat*10+1)
    hold on
    plot(paretoSats{:},paretoPlanes{:},'--*')
    xlabel('# Sats')
    ylabel('# Planes')
    xlim([50,80])
    ylim([0,60])
    grid on
    legend('Optimal','Nominal')
    hold off
end