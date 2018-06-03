function [ fit ] = WalkerFitnessOscPat(x, nSatsT, nPlanesP, phasingF, timeVec,...
                                       latGS, lonGS, elevMin, relTol, absTol)
%WalkerFitnessOscPat Simulates a Walker Constellation and calculates the
%fitness.  Used for patternsearch where T,P,F are given

% Initialization
inc = x(1);
alt = x(2);
gmst0 = x(3);

WC = WalkerConstellation(nSatsT,nPlanesP,phasingF,inc,alt);
Prop = Propagator(WC,relTol,absTol);

% Propagate
[propTime,xEci] = Prop.PropEciJ2(timeVec);

% Evaluate Performance
pdop = TdoaPdopVec(xEci,propTime,latGS,lonGS,gmst0,elevMin);
fit = max(pdop);
end