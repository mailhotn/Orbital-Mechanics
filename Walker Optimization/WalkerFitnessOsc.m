function [ fit ] = WalkerFitnessOsc(x, nSatsT, timeVec, latGs, lonGs,...
                                    elevMin, relTol, absTol)
%WalkerFitnessOsc Simulates a Walker Constellation and calculates the
%fitness

% Initialization
nPlanesP = x(1);
phasingF = x(2);
inc      = x(3);
alt      = x(4);
gmst0 = x(5);

WC = WalkerConstellation(nSatsT,nPlanesP,phasingF,inc,alt);
Prop = Propagator(WC,relTol,absTol);

% Propagate
[simTime,eciState] = Prop.PropEciJ2(timeVec);

% Evaluate Performance
pdop = TdoaPdopVec(eciState,simTime,latGs,lonGs,gmst0,elevMin);
fit = max(pdop);
end