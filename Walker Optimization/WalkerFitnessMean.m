function [ fit ] = WalkerFitnessMean(x, nSatsT, timeVec,...
                                     latGs, lonGs, elevMin, relTol, absTol)
%WalkerFitnessMean Simulates a Walker Constellation and calculates the
%fitness

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
xEci  = reshape([R;V],6*nSatsT,length(propTime));
xEci  = xEci.';
% Evaluate Performance

pdop = TdoaPdopVec(xEci,propTime,latGs,lonGs,gmst0,elevMin);
fit = max(pdop);
end