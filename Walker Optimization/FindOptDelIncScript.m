datafolder = 'C:\Users\User\Dropbox\Walker Optimization Data\Previous Optimization Runs\Walker lat30-40 dInc test';
%% Set Performance Goals
intTarget = 10;
maxTarget = 1000;
covTarget = 99;
p90Target = 5;
load([datafolder '\OptParams.mat']);

nSats = minSats:maxSats;
maxPdop = nan(numel(hAList),numel(nSats));
intPdop = nan(numel(hAList),numel(nSats));
coverage = nan(numel(hAList),numel(nSats));
p90 = nan(numel(hAList),numel(nSats));
for iLat = 1:numel(latList)
    paretoSats = inf;
    paretoPlanes = inf;
    paretoDelInc = inf;
    for iSats = 1:numel(nSats)
        for iHA = numel(hAList):-1:1
            load([datafolder '\LatticeExSol_Lat_'...
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
                        paretoSats = [paretoSats, nSats(iSats)];
                        paretoPlanes = [paretoPlanes, nPlanesToAchieve];
                        paretoDelInc = [paretoDelInc, hAList(iHA)];
                    else
                        paretoSats = [nSats(iSats)];
                        paretoPlanes = [nPlanesToAchieve];
                        paretoDelInc = [hAList(iHA)];
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
    plot(paretoSats,paretoPlanes,'--o')
    text(paretoSats,paretoPlanes,num2str(latList(iLat)+paretoDelInc(:)))
end