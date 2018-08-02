function [ fit ] = FlowerFitnessMean( x, nSats, PropParams)
%FlowerFitnessMean props an FC and calculates the fitness for GA
%
% Propagation is done with the mean elements. Fitness is calculated based
% on the PDOP relative to a given ground station.
% Genome(9): nPetals-fN-fD-fH-w-i-altP-raan0-M0
%
% PropParams = {latGs, lonGs, elevMin, relTol, absTol, maxPdop, fitFunc}
%% Repack Genome & Initialize Constellation
Architecture.nPetals = x(1);
Architecture.nSats   = nSats;

Phasing.fN     = x(2);
Phasing.fD     = x(3);
Phasing.fH     = x(4);

Orbits.w        = x(5);
Orbits.inc      = x(6);
Orbits.altP     = x(7);

Initial.raan0  = x(8);
Initial.M0     = x(9);

Architecture.nDays   = Architecture.nSats/Phasing.fD;

FC = FlowerConstellation(Architecture,Phasing,Orbits,Initial);
Prop = Propagator(FC, PropParams.relTol, PropParams.absTol);
timeVec = 0:100:86164*FC.nDays;
%% Propagate
[propTime,propState] = Prop.PropOeMeanFast(timeVec);

% Transform to ECI
oeMean = reshape(propState.',6,length(propTime)*nSats);
oeMean(6,:) = me2ta(oeMean(6,:),oeMean(2,:));
[R, V] = oe2eci(oeMean);
xEci  = reshape([R;V],6*nSats,length(propTime));
xEci  = xEci.';

% Evaluate Performance
[pdop, ~] = TdoaPdopVec(xEci,propTime,PropParams.latGs,PropParams.lonGs,...
                        0,PropParams.elevMin);
switch PropParams.fitFunc
    case 'Integral'
        pdop(pdop > PropParams.maxPdop) = PropParams.maxPdop;
        pdop(isnan(pdop)) = PropParams.maxPdop;
        fit = trapz(propTime,pdop)/(propTime(end)-propTime(1));
    case 'Mean'
        pdop(pdop > PropParams.maxPdop) = PropParams.maxPdop;
        pdop(isnan(pdop)) = PropParams.maxPdop;
        fit = mean(pdop);
    case 'Max'
        pdop(pdop > PropParams.maxPdop) = PropParams.maxPdop;
        pdop(isnan(pdop)) = PropParams.maxPdop*10;
        fit = max(pdop);
    otherwise
        fit = inf;
end