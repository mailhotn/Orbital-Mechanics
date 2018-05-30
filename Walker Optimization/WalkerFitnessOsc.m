function [ fit ] = WalkerFitnessOsc(x, nSatsT, timeVec, latGs, lonGs,...
    elevMin, relTol, absTol)
%WalkerFitness_sphere Simulates a Walker Constellation and calculates the
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
[simTime,eciState] = Prop.prop_ECI_J2(timeVec);

% Evaluate Performance
pdop = get_PDOP_vec_WGS84(eciState,simTime,latGs,lonGs,gmst0,elevMin);
fit = max(pdop);
end