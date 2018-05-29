function [ fit ] = WalkerFitnessOscPat(x, T, P, F, time, latGS, lonGS, eMin, relTol, absTol)
%WalkerFitness_sphere Simulates a Walker Constellation and calculates the
%fitness

% Initialization
inclination = x(1);
altitude = x(2);
GMST0 = x(3);

WC = WalkerConstellation(T,P,F,inclination,altitude);
Prop = Propagator(WC,relTol,absTol);

% Propagate
[tVec,X_ECI] = Prop.prop_ECI_J2(time);

% Evaluate Performance
PDOP = get_PDOP_vec_WGS84(X_ECI,tVec,latGS,lonGS,GMST0,eMin);
fit = max(PDOP);
end