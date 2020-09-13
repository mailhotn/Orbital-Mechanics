function [fitness] = LatticeGaFitness(x, Arch, latEm, OptParams)

% Handle Input
hA = OptParams.hAList(x(1));
inc = latEm + x(2);
Phase = struct();
Phase.nC1 = x(3);
Phase.nC2 = x(4);
Phase.nC3 = x(5);

% Orbital Parameters
[sma,ecc] = CalcRgtSmaApoHeight(inc, hA, OptParams.nRepeats, OptParams.nDays);
Orbit = struct();
Orbit.sma = sma;
Orbit.ecc = ecc;
Orbit.inc = inc;
Orbit.hA = hA;

% Initial Conditions
InitCon = InitConElliptical(ecc,inc,sma,latEm,0);

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

fitness = (2*intPdopNom + intPdopNorth + intPdopSouth)/4 + ...
    abs(intPdopNorth - intPdopSouth);


