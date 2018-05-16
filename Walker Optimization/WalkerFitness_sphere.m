function [ fit ] = WalkerFitness_sphere(x,T)
%WalkerFitness_sphere Simulates a Walker Constellation and calculates the
%fitness
%   Simulates a constellation for several days then calculates the PDOP every
%   100 seconds. The fitness is then taken calculated from the PDOP.

% Initialization
P = x(1);
F = x(2);
inclination = x(3);
altitude = x(4);
GMST0 = x(5);

WC = WalkerConstellation(T,P,F,inclination,altitude);
Prop = Propagator(WC,1e-6,1e-6);
Time = 0:100:5*86400;

% Propagate
[T_vec,X] = Prop.prop_ECI_J2(Time);

% Evaluate Performance
PDOP = get_PDOP_vec_sphere(X,T_vec,30,30,GMST0,15);

fit = (PDOP.'*PDOP)/length(PDOP) + max(PDOP)^2;
end