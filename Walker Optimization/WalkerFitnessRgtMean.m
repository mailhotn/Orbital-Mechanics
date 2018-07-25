function [ fit ] = WalkerFitnessRgtMean(x, nSatsT, timeVec,latGs, lonGs,...
                                        elevMin, relTol, absTol,maxPdop,...
                                        fitFunc)
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
oeMean = reshape(propState.',6,length(propTime)*nSatsT);
oeMean(6,:) = me2ta(oeMean(6,:),oeMean(2,:));
[R, V] = oe2eci(oeMean);
xEci  = reshape([R;V],6*nSatsT,length(propTime));
xEci  = xEci.';

% Evaluate Performance
[pdop, ~] = TdoaPdopVec(xEci,propTime,latGs,lonGs,gmst0,elevMin);
switch fitFunc
    case 'Integral'
        pdop(pdop > maxPdop) = maxPdop;
        pdop(isnan(pdop)) = maxPdop;
        fit = trapz(propTime,pdop)/(propTime(end)-propTime(1));
    case 'Mean'
        pdop(pdop > maxPdop) = maxPdop;
        pdop(isnan(pdop)) = maxPdop;
        fit = mean(pdop);
    case 'Max'
        pdop(pdop > maxPdop) = maxPdop;
        pdop(isnan(pdop)) = maxPdop;
        fit = max(pdop);
    otherwise
        fit = inf;
end