function [fitness] = LatMoGaFitness(x, latEm, OptParams)

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
% Remove Jumps
pdop(pdop > OptParams.maxPdop) = OptParams.maxPdop;
pdop(isnan(pdop)) = OptParams.maxPdop;
pdopN(pdopN > OptParams.maxPdop) = OptParams.maxPdop;
pdopN(isnan(pdopN)) = OptParams.maxPdop;
pdopS(pdopS > OptParams.maxPdop) = OptParams.maxPdop;
pdopS(isnan(pdopS)) = OptParams.maxPdop;

% Integrate
intPdopNom = trapz(propTime,pdop)/(propTime(end)-propTime(1));
intPdopNorth = trapz(propTime,pdopN)/(propTime(end)-propTime(1));
intPdopSouth = trapz(propTime,pdopS)/(propTime(end)-propTime(1));

% fitness
fitness = inf(1,3);
fitness(1) = x(1);
fitness(2) = x(2);
fitness(3) = min([(intPdopNom + intPdopNorth + intPdopSouth)/3,20]);
