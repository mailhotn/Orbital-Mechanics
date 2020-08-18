function [ ExSol ] = LatticeExSearch( Arch, Orbit, InitCon, latGs, PropParams)
%LatticeExSearch Exhaustively checks nPlanes, nC1, nC2, nC3 for the best ones
% Constrained to have a minimal number of seperate ground tracks.
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
            Arch.nSatsPerAop = Arch.nDays; % Constraint!  can be different if nDays > 1
            Arch.nAops = Arch.nSats/Arch.nPlanes/Arch.nSatsPerAop;
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
                        minDist = CalcMinDist(LC);
                        if minDist > PropParams.minMinDist
                            
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