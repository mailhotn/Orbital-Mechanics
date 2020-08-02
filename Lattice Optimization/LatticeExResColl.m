function [ExSol,nColls] = LatticeExResColl(latEm,nSats,hA,PropParams,...
    sourceFolder,targetFolder)
%LatticeExResColl Fixes constellations with collisions
%   checks existing solutions for collision and redoes optimization if any
%   are found.
nColls = 0;
load([sourceFolder '\LatticeExSol_Lat_' num2str(latEm) '_nSats_' ...
    num2str(nSats) '_hA_' num2str(hA) '.mat'],'ExSol');

minDists = zeros(1,size(ExSol.archMat,2));

for iPlanes = 1:size(ExSol.archMat,2)
    Arch.nSats = nSats;
    Arch.nPlanes = ExSol.archMat(1,iPlanes);
    Arch.nAops = ExSol.archMat(2,iPlanes);
    Arch.nSatsPerAop = ExSol.archMat(3,iPlanes);
    Arch.nRepeats = 14;
    Arch.nDays = 1;
    
    Phase.nC1 = ExSol.phaseMat(1,iPlanes);
    Phase.nC2 = ExSol.phaseMat(2,iPlanes);
    Phase.nC3 = ExSol.phaseMat(3,iPlanes);
    
    Orbit = ExSol.orbits{iPlanes};
    InitCon = ExSol.inits{iPlanes};
    
    LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);
    
    minDist = CalcMinDist(LC);
    
    if minDist < PropParams.minMinDist
        nColls = nColls + 1; % Collision Detected!!!
        
        % Redo optimization with collision checking
        
        incList = -1:0.2:5;
        ExSol.intPdop(iPlanes) = inf;
        
        for iInc = 1:length(incList)
            inc = latEm + 10 + incList(iInc);
            [sma, ecc] = CalcRgtSmaApoHeight(inc,hA,Arch.nRepeats, Arch.nDays);
            Orbit.sma = sma;
            Orbit.ecc = ecc;
            Orbit.inc = inc;
            
            InitCon = InitConElliptical(Orbit.ecc,Orbit.inc,Orbit.sma,latEm,0);
            
            nC3List = divisors(Arch.nPlanes); % Constraint! nC3 does not have to be a divisor of nPlanes
            for iNC3 = 1:length(nC3List)
                Phase.nC3 = nC3List(iNC3);
                for nC2 = 1:Arch.nAops
                    Phase.nC2 = nC2;
                    gcdO = gcd(Arch.nPlanes,Phase.nC3);
                    lamList = ceil((Arch.nRepeats-Arch.nPlanes)/gcdO):1:...
                        floor((Arch.nRepeats-1)/gcdO);
                    for iLam = 1:length(lamList)
                        Phase.nC1 = Arch.nRepeats - lamList(iLam)*gcdO;
                        % Create Constellation and propagate
                        LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);
                        
                        minDist = CalcMinDist(LC); %Collision check
                        if minDist > PropParams.minMinDist
                            Prop = Propagator(LC,PropParams.relTol,PropParams.absTol);
                            [propTime, propState] = Prop.PropEciJ2(PropParams.timeVec);
                            % Evaluate PDOP
                            [pdopC, ~] = TdoaPdopVec(propState,propTime,latEm,0,0,PropParams.elevMin);
                            [pdopN, ~] = TdoaPdopVec(propState,propTime,...
                                latEm+PropParams.delLat,0,0,PropParams.elevMin);
                            [pdopS, ~] = TdoaPdopVec(propState,propTime,...
                                latEm-PropParams.delLat,0,0,PropParams.elevMin);
                            pdop = [pdopN.';pdopC.';pdopS.'];
                            % Calculate Performance criteria
                            coverage = 100 - sum(isnan(pdop),2)/length(pdop)*100;
                            if any(~isnan(pdop),'all')
                                maxPdop  = max(pdop(~isnan(pdop)),[],2);
                                pdop(pdop > PropParams.maxPdop) = PropParams.maxPdop;
                                pdop(isnan(pdop)) = PropParams.maxPdop;
                                intPdop  = trapz(propTime,pdop,2)/(propTime(end)-propTime(1));
                                intPdop = mean(intPdop,1);
                                if intPdop < ExSol.intPdop(iPlanes)
                                    ExSol.intPdop(iPlanes) = intPdop;
                                    ExSol.maxPdop(iPlanes) = mean(maxPdop,1);
                                    ExSol.coverage(iPlanes) = mean(coverage,1);
                                    ExSol.p90(iPlanes) = mean(prctile(pdop,90,2),1);
                                    ExSol.p75(iPlanes) = mean(prctile(pdop,75,2),1);
                                    ExSol.p50(iPlanes) = mean(prctile(pdop,50,2),1);
                                    ExSol.phaseMat(:,iPlanes) = [Phase.nC1;Phase.nC2;Phase.nC3];
                                    ExSol.archMat(:,iPlanes)  = [Arch.nPlanes;Arch.nAops;Arch.nSatsPerAop];
                                    ExSol.orbits{iPlanes} = Orbit;
                                    ExSol.inits{iPlanes} = InitCon;
                                end
                            end
                        end
                    end
                end
            end
        end
    else
        minDists(iPlanes) = minDist;
    end
end
[fit,iOpt] = min(ExSol.intPdop);
ExSol.iOpt = iOpt;
ExSol.fit = fit;
ExSol.minDists = minDists;
save([targetFolder '\LatticeExSol_Lat_' num2str(latEm)...
        '_nSats_' num2str(nSats) '_hA_' num2str(hA) '.mat'],'ExSol');


