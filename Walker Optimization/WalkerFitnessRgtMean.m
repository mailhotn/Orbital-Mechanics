function [ fit ] = WalkerFitnessRgtMean(x, nSatsT, timeVec,...
                                     latGs, lonGs, elevMin, relTol, absTol)
%WalkerFitnessRgtMean Simulates a Walker Constellation and calculates the
%fitness

% Initialization
nPlanesP  = x(1);
phasingF  = x(2);
nRepeatsj = x(3);
inc = x(4);
gmst0 = x(5);
sma = CalcRgtElement([],0,inc,nRepeatsj,1);

WC = WalkerConstellation(nSatsT, nPlanesP, phasingF, inc, sma - 6378.137);
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

[pdop, ~] = TdoaPdopVec(xEci,propTime,latGs,lonGs,gmst0,elevMin);

pdop(isnan(pdop)) = max(pdop)*10;
fit = trapz(propTime,pdop)/(propTime(end)-propTime(1));
end