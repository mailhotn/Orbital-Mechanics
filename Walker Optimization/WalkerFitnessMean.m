function [ fit ] = WalkerFitnessMean(x, nSatsT, timeVec,...
                                     latGs, lonGs, elevMin, relTol, absTol)
%WalkerFitness_sphere Simulates a Walker Constellation and calculates the
%fitness
%   Simulates a constellation for several days then calculates the PDOP every
%   100 seconds. The fitness is then taken calculated from the PDOP.

% Initialization
nPlanesP = x(1);
phasingF = x(2);
inc = x(3);
alt = x(4);
gmst0 = x(5);

WC = WalkerConstellation(nSatsT, nPlanesP, phasingF, inc, alt);
Prop = Propagator(WC, relTol, absTol);

% Propagate
[propTime,propState] = Prop.PropOeMeanFast(timeVec);

% Transform to ECI
oeOsc = me2osc(reshape(propState.',6,length(propTime)*nSatsT));
oeOsc(6,:) = me2ta(oeOsc(6,:),oeOsc(2,:));
[R, V] = oe2eci(oeOsc);
X_ECI  = reshape([R;V],6*nSatsT,length(propTime));
X_ECI  = X_ECI.';
% Evaluate Performance

PDOP = get_PDOP_vec_WGS84(X_ECI,propTime,latGs,lonGs,gmst0,elevMin);
fit = max(PDOP);
end