function [ ExSol ] = LatticeExSearch( Arch, Orbit, InitCon, latGs, PropParams)
%LatticeExSearch Exhaustively checks nPlanes, nC2, for the best ones
%   Output Matrices columns correspond to different values of nPlanes
%   Rows correspond to different values of nC2
%   nPlanes --->
%       nC2 || x x x x 
%           || x x x x 
%           \/ x x x x 

% initialize arrays
pList = divisors(Arch.nSats);
meanMat = nan(max(pList));
maxMat  = nan(max(pList));
intMat  = nan(max(pList));
covMat  = nan(max(pList));

for iPlanes = 1:length(pList)
    Arch.nPlanes = pList(iPlanes);
    Phase.nC3 = Arch.nPlanes; % Constraint! can be different, would increase number of relative orbits.
    Arch.nSatsPerAop = Arch.nDays; % Constraint!  can be different if nDays > 1
    Arch.nAops = Arch.nSats/Arch.nPlanes/Arch.nSatsPerAop;
    gcdOrbits = gcd(Arch.nPlanes,Phase.nC3);
    if gcdOrbits > Arch.nRepeats
        Phase.nC1 = Arch.nRepeats;
    else
        l = floor(Arch.nRepeats/gcdOrbits);
        if l*gcdOrbits == Arch.nRepeats
            l = l-1;
        end
        Phase.nC1 = Arch.nRepeats - l*gcdOrbits;
    end
    for nC2 = 1:Arch.nAops
        Phase.nC2 = nC2;
        % Create Constellation and propagate
        LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);
        Prop = Propagator(LC,PropParams.relTol,PropParams.absTol);
        [propTime, propState] = Prop.PropEciJ2(PropParams.timeVec);
        % Evaluate PDOP
        [pdop, ~] = TdoaPdopVec(propState,propTime,latGs,0,0,PropParams.elevMin);
        % Calculate Performance criteria
        covMat(nC2,pList(iPlanes))  = 100 - sum(isnan(pdop))/length(pdop)*100;
        if any(~isnan(pdop))
            maxMat(nC2,pList(iPlanes))  = max(pdop(~isnan(pdop)));
            pdop(pdop > PropParams.maxPdop)    = PropParams.maxPdop;
            pdop(isnan(pdop)) = PropParams.maxPdop;
            meanMat(nC2,pList(iPlanes)) = mean(pdop);
            intMat(nC2,pList(iPlanes))  = trapz(propTime,pdop)/...
                (propTime(end)-propTime(1));
        end
    end
end
% Find Best Solution
[minForP,indP] = min(intMat);
[~,optNPlanes] = min(minForP);
optNC2 = indP(optNPlanes);
% Assign Values to Output Struct
ExSol.Arch  = Arch;
ExSol.Phase = Phase;
ExSol.Orbit = Orbit;
ExSol.InitCon = InitCon;
ExSol.latGs = latGs;
ExSol.optNPlanes  = optNPlanes;
ExSol.optNC2  = optNC2;
ExSol.PropParams = PropParams;
ExSol.coverage = covMat;
ExSol.maxPdop  = maxMat;
ExSol.meanPdop = meanMat;
ExSol.intPdop  = intMat;
ExSol.fit = ExSol.intPdop(ExSol.optNC2, ExSol.optNPlanes);
save([PropParams.datafolder '\LatticeExSol_Lat_' num2str(latGs)...
                '_nSats_' num2str(Arch.nSats) '_ecc_' num2str(Orbit.ecc) '.mat'],'ExSol');
end