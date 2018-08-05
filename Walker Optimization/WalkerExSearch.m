function [ ExSol ] = WalkerExSearch( nSats, inc, sma, raan0, latGs, PropParams)
%WalkerExSearch Exhaustively checks P & F for the best ones
%   Detailed explanation goes here
primary = earth();
% initialize arrays
pList = divisors(nSats);
meanMat = nan(max(pList));
maxMat  = nan(max(pList));
intMat  = nan(max(pList));
covMat  = nan(max(pList));
% Go over all P anf F
for iPlanes = 1:length(pList)
    for phasingF = 0:(pList(iPlanes)-1)
        % Create Constellation and propagate
        WC = WalkerConstellation(nSats, pList(iPlanes), phasingF, inc, ...
                                 sma-primary.Re,raan0,primary);
        Prop = Propagator(WC,PropParams.relTol,PropParams.absTol);
        [propTime, propState] = Prop.PropEciJ2(PropParams.timeVec);
        % Evaluate PDOP
        [pdop, ~] = TdoaPdopVec(propState,propTime,latGs,0,0,PropParams.elevMin);
        % Calculate Performance criteria
        covMat(phasingF+1,pList(iPlanes))  = 100 - sum(isnan(pdop))/length(pdop)*100;
        if any(~isnan(pdop))
            maxMat(phasingF+1,pList(iPlanes))  = max(pdop(~isnan(pdop)));
            pdop(pdop > PropParams.maxPdop)    = PropParams.maxPdop;
            pdop(isnan(pdop)) = PropParams.maxPdop;
            meanMat(phasingF+1,pList(iPlanes)) = mean(pdop);
            intMat(phasingF+1,pList(iPlanes))  = trapz(propTime,pdop)/...
                (propTime(end)-propTime(1));
        end
    end
end
% Find Best Solution
[minForP,indP] = min(intMat);
[~,optP] = min(minForP);
optF = indP(optP) - 1;
% Assign Values to Output Struct
ExSol.nSats = nSats;
ExSol.latGs = latGs;
ExSol.optP  = optP;
ExSol.optF  = optF;
ExSol.inc   = inc;
ExSol.alt   = sma-primary.Re;
ExSol.raan0 = raan0;
ExSol.PropParams = PropParams;
ExSol.coverage = covMat;
ExSol.maxPdop = maxMat;
ExSol.meanPdop = meanMat;
ExSol.intPdop = intMat;
ExSol.fit = ExSol.intPdop(ExSol.optF + 1,ExSol.optP);

save([PropParams.datafolder '\WalkerRgtExSol_Lat_' num2str(latGs)...
                '_T_' num2str(nSats) '.mat'],'ExSol');
end