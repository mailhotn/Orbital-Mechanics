function [ ExSol ] = LatticeExSearch( Arch, Orbit, InitCon, latGs, PropParams)
%LatticeExSearch Exhaustively checks nPlanes, nC1, nC2, nC3 for the best ones
% Constrained to have a minimal number of seperate ground tracks.
pList = divisors(Arch.nSats);

% initialize Performance Arrays
maxVec = inf(1,length(pList));
intVec = inf(1,length(pList));
covVec = zeros(1,length(pList));

% Initialize Constellation Parameter Arrays
phaseParams = nan(3,length(pList)); % nC1, nC2, nC3
archParams  = nan(3,length(pList)); % nPlanes, nAops, nSatsPerAop

for iPlanes = 1:length(pList)
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
                Prop = Propagator(LC,PropParams.relTol,PropParams.absTol);
                [propTime, propState] = Prop.PropEciJ2(PropParams.timeVec);
                % Evaluate PDOP
                [pdop, ~] = TdoaPdopVec(propState,propTime,latGs,0,0,PropParams.elevMin);
                % Calculate Performance criteria
                coverage = 100 - sum(isnan(pdop))/length(pdop)*100;
                if any(~isnan(pdop))
                    maxPdop  = max(pdop(~isnan(pdop)));
                    pdop(pdop > PropParams.maxPdop)    = PropParams.maxPdop;
                    pdop(isnan(pdop)) = PropParams.maxPdop;
                    intPdop  = trapz(propTime,pdop)/(propTime(end)-propTime(1));
                    if intPdop < intVec(iPlanes)
                        intVec(iPlanes) = intPdop;
                        maxVec(iPlanes) = maxPdop;
                        covVec(iPlanes) = coverage;
                        phaseParams(:,iPlanes) = [Phase.nC1;Phase.nC2;Phase.nC3];
                        archParams(:,iPlanes)  = [Arch.nPlanes;Arch.nAops;Arch.nSatsPerAop];
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
ExSol.Orbit = Orbit;
ExSol.InitCon = InitCon;
ExSol.latGs = latGs;
ExSol.optNPlanes = archParams(1,iOpt);
ExSol.iOpt = iOpt;
ExSol.PropParams = PropParams;
ExSol.coverage = covVec;
ExSol.maxPdop  = maxVec;
ExSol.intPdop  = intVec;
ExSol.fit = fit;
save([PropParams.datafolder '\LatticeExSol_Lat_' num2str(latGs)...
    '_nSats_' num2str(Arch.nSats) '_ecc_' num2str(Orbit.ecc) '.mat'],'ExSol');
end