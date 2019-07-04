function [GaSol] = LatMoGaFitnessPlus(x, latEm, OptParams)

% Handle input
Arch = struct();
Arch.nSats = x(1);
Arch.nPlanes = x(2);
Arch.nSatsPerAop = 1;
Arch.nAops = Arch.nSats/Arch.nPlanes/Arch.nSatsPerAop;

Phase = struct();
Phase.nC1 = x(3);
Phase.nC2 = x(4);
Phase.nC3 = x(5);

[sma,ecc] = CalcRgtSmaApoHeight(x(6), x(7), OptParams.nRepeats, OptParams.nDays);
Orbit = struct();
Orbit.sma = sma;
Orbit.ecc = ecc;
Orbit.inc = x(6);
Orbit.hA = x(7);

InitCon = InitConElliptical(Orbit.ecc,Orbit.inc,Orbit.sma,latEm,0);
% initialize Performance Arrays
maxMat = inf(3,1);
intMat = inf(3,1);
p90Mat = inf(3,1);
p75Mat = inf(3,1);
p50Mat = inf(3,1);
covMat = zeros(3,1);

% Creat Constellation & Propagate
LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);
Prop = Propagator(LC,OptParams.relTol,OptParams.absTol);
[propTime, propState] = Prop.PropEciJ2(OptParams.timeVec);

% Evaluate PDOP
[pdop, ~] = TdoaPdopVec(propState,propTime,latEm,0,0,OptParams.elevMin);
[pdopN, ~] = TdoaPdopVec(propState,propTime,...
    latEm+OptParams.delLat,0,0,OptParams.elevMin);
[pdopS, ~] = TdoaPdopVec(propState,propTime,...
    latEm-OptParams.delLat,0,0,OptParams.elevMin);

pdopMat = [pdopN.';pdop.';pdopS.'];
% Evaluate Performance at all latitudes
for iLat = 1:3
    coverage = 100 - sum(isnan(pdopMat(iLat,:)))/length(pdopMat(iLat,:))*100;
    if any(~isnan(pdopMat(iLat,:)))
        maxPdop  = max(pdopMat(iLat,~isnan(pdopMat(iLat,:))));
        pdopMat(iLat,pdopMat(iLat,:) > OptParams.maxPdop) = OptParams.maxPdop;
        pdopMat(iLat,isnan(pdopMat(iLat,:))) = OptParams.maxPdop;
        intPdop  = trapz(propTime,pdopMat(iLat,:))/(propTime(end)-propTime(1));
        intMat(iLat) = intPdop;
        maxMat(iLat) = maxPdop;
        covMat(iLat) = coverage;
        p90Mat(iLat) = prctile(pdopMat(iLat,:),90);
        p75Mat(iLat) = prctile(pdopMat(iLat,:),75);
        p50Mat(iLat) = prctile(pdopMat(iLat,:),50);
    end
end

GaSol = struct();
GaSol.Con = LC;
GaSol.coverage = covMat;
GaSol.maxPdop  = maxMat;
GaSol.intPdop  = intMat;
GaSol.p90 = p90Mat;
GaSol.p75 = p75Mat;
GaSol.p50 = p50Mat;
