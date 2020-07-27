function [ ExSol ] = WalkerExSearch( Arch, Orbit, InitCon, latGs, PropParams)
%WalkerExSearch Exhaustively checks P & F for the best ones
%   Actually Lattice Flower Constellations
if ~exist([PropParams.datafolder '\LatticeExSol_Lat_' num2str(latGs)...
        '_nSats_' num2str(Arch.nSats) '_hA_' num2str(Orbit.hA) '.mat'],'file')
    pList = divisors(Arch.nSats);
    incList = -1:0.2:5;
    Orbit0 = Orbit;
    
    % initialize Performance Arrays
    maxVec = inf(1,length(pList));
    intVec = inf(1,length(pList));
    p90Vec = inf(1,length(pList));
    p75Vec = inf(1,length(pList));
    p50Vec = inf(1,length(pList));
    covVec = zeros(1,length(pList));
    
    % Initialize Constellation Parameter Arrays
    phaseParams = nan(3,length(pList)); % nC1, nC2, nC3
    archParams  = nan(3,length(pList)); % nPlanes, nAops, nSatsPerAop
    Phase.nC3 = 0;
    Phase.nC2 = 0;
    orbits = cell(1,length(pList));
    inits  = cell(1,length(pList));
    
    for iPlanes = 1:length(pList)
        for iInc = 1:length(incList)
            % Modify Inclination & redo initcon
            inc = Orbit0.inc + incList(iInc);
            [sma, ecc] = CalcRgtSmaApoHeight(inc,Orbit.hA,Arch.nRepeats, Arch.nDays);
            Orbit.sma = sma;
            Orbit.ecc = ecc;
            Orbit.inc = inc;
            
            InitCon = InitConElliptical(Orbit.ecc,Orbit.inc,Orbit.sma,latGs,0);
            
            Arch.nPlanes = pList(iPlanes);
            Arch.nSatsPerAop = Arch.nSats/Arch.nPlanes;
            Arch.nAops = 1;
            for nC1 = 0:Arch.nPlanes-1
                Phase.nC1 = -nC1;
                % Create Constellation and propagate
                LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);
                Prop = Propagator(LC,PropParams.relTol,PropParams.absTol);
                [propTime, propState] = Prop.PropEciJ2(PropParams.timeVec);
                % Evaluate PDOP
                [pdopC, ~] = TdoaPdopVec(propState,propTime,latGs,0,0,PropParams.elevMin);
                [pdopN, ~] = TdoaPdopVec(propState,propTime,...
                    latGs+PropParams.delLat,0,0,PropParams.elevMin);
                [pdopS, ~] = TdoaPdopVec(propState,propTime,...
                    latGs-PropParams.delLat,0,0,PropParams.elevMin);
                pdop = [pdopN.';pdopC.';pdopS.'];
                % Calculate Performance criteria
                coverage = 100 - sum(isnan(pdop),2)/length(pdop)*100;
                if any(~isnan(pdop),'all')
                    maxPdop  = max(pdop(~isnan(pdop)),[],2);
                    pdop(pdop > PropParams.maxPdop) = PropParams.maxPdop;
                    pdop(isnan(pdop)) = PropParams.maxPdop;
                    intPdop  = trapz(propTime,pdop,2)/(propTime(end)-propTime(1));
                    intPdop = mean(intPdop,1);
                    if intPdop < intVec(iPlanes)
                        intVec(iPlanes) = intPdop;
                        maxVec(iPlanes) = mean(maxPdop,1);
                        covVec(iPlanes) = mean(coverage,1);
                        p90Vec(iPlanes) = mean(prctile(pdop,90,2),1);
                        p75Vec(iPlanes) = mean(prctile(pdop,75,2),1);
                        p50Vec(iPlanes) = mean(prctile(pdop,50,2),1);
                        phaseParams(:,iPlanes) = [Phase.nC1;Phase.nC2;Phase.nC3];
                        archParams(:,iPlanes)  = [Arch.nPlanes;Arch.nAops;Arch.nSatsPerAop];
                        orbits{iPlanes} = Orbit;
                        inits{iPlanes} = InitCon;
                    end
                end
            end
        end
    end
    % Find Best Solution
    [fit,iOpt] = min(intVec);
    % Assign Values to Output Struct
    ExSol.nSats = Arch.nSats;
    ExSol.archMat  = archParams;
    ExSol.phaseMat = phaseParams;
    ExSol.orbits = orbits;
    ExSol.inits = inits;
    ExSol.latGs = latGs;
    ExSol.optNPlanes = archParams(1,iOpt);
    ExSol.iOpt = iOpt;
    ExSol.PropParams = PropParams;
    ExSol.coverage = covVec;
    ExSol.maxPdop  = maxVec;
    ExSol.intPdop  = intVec;
    ExSol.p90 = p90Vec;
    ExSol.p75 = p75Vec;
    ExSol.p50 = p50Vec;
    ExSol.fit = fit;
    save([PropParams.datafolder '\LatticeExSol_Lat_' num2str(latGs)...
        '_nSats_' num2str(Arch.nSats) '_hA_' num2str(Orbit.hA) '.mat'],'ExSol');
else
    load([PropParams.datafolder '\LatticeExSol_Lat_' num2str(latGs)...
        '_nSats_' num2str(Arch.nSats) '_hA_' num2str(Orbit.hA) '.mat']);
end


%% Old non-lattice Version
% primary = earth();
% % initialize arrays
% pList = divisors(nSats);
% meanMat = nan(max(pList));
% maxMat  = nan(max(pList));
% intMat  = nan(max(pList));
% covMat  = nan(max(pList));
% p90Mat  = nan(max(pList));
% p75Mat  = nan(max(pList));
% p50Mat  = nan(max(pList));
% % Go over all P anf F
% for iPlanes = 1:length(pList)
%     for phasingF = 0:(pList(iPlanes)-1)
%         % Create Constellation and propagate
%         WC = WalkerConstellation(nSats, pList(iPlanes), phasingF, inc, ...
%                                  sma-primary.Re,raan0,primary);
%         Prop = Propagator(WC,PropParams.relTol,PropParams.absTol);
%         [propTime, propState] = Prop.PropEciJ2(PropParams.timeVec);
%         % Evaluate PDOP
%         [pdop, ~] = TdoaPdopVec(propState,propTime,latGs,0,0,PropParams.elevMin);
%         % Calculate Performance criteria
%         covMat(phasingF+1,pList(iPlanes))  = 100 - sum(isnan(pdop))/length(pdop)*100;
%         if any(~isnan(pdop))
%             maxMat(phasingF+1,pList(iPlanes))  = max(pdop(~isnan(pdop)));
%             pdop(pdop > PropParams.maxPdop)    = PropParams.maxPdop;
%             pdop(isnan(pdop)) = PropParams.maxPdop;
%             meanMat(phasingF+1,pList(iPlanes)) = mean(pdop);
%             intMat(phasingF+1,pList(iPlanes))  = trapz(propTime,pdop)/...
%                 (propTime(end)-propTime(1));
%             p90Mat(phasingF+1,pList(iPlanes)) = prctile(pdop,90);
%             p75Mat(phasingF+1,pList(iPlanes)) = prctile(pdop,75);
%             p50Mat(phasingF+1,pList(iPlanes)) = prctile(pdop,50);
%         end
%     end
% end
% % Find Best Solution
% [minForP,indP] = min(intMat);
% [~,optP] = min(minForP);
% optF = indP(optP) - 1;
% % Assign Values to Output Struct
% ExSol.nSats = nSats;
% ExSol.latGs = latGs;
% ExSol.optP  = optP;
% ExSol.optF  = optF;
% ExSol.inc   = inc;
% ExSol.alt   = sma-primary.Re;
% ExSol.raan0 = raan0;
% ExSol.PropParams = PropParams;
% ExSol.coverage = covMat;
% ExSol.maxPdop = maxMat;
% ExSol.meanPdop = meanMat;
% ExSol.intPdop = intMat;
% ExSol.p90 = p90Mat;
% ExSol.p75 = p75Mat;
% ExSol.p50 = p50Mat;
% ExSol.fit = ExSol.intPdop(ExSol.optF + 1,ExSol.optP);
%
% save([PropParams.datafolder '\WalkerRgtExSol_Lat_' num2str(latGs)...
%                 '_T_' num2str(nSats) '.mat'],'ExSol');
% end