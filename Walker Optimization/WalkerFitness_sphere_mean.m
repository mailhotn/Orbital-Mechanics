function [ fit ] = WalkerFitness_sphere_mean(x, T, Time)
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
% Time = 0:1000:30*86400;

% Propagate

[T_vec,X] = Prop.prop_OE_Mean_lin(Time);

% Transform to ECI
tic
OE_osc = me2osc(reshape(X.',6,length(T_vec)*T));
toc
OE_osc(6,:) = me2ta(OE_osc(6,:),OE_osc(2,:));
tic
[R, V] = oe2eci(OE_osc);
toc
X_ECI  = reshape([R;V],6*T,length(T_vec));
X_ECI  = X_ECI.';
% Evaluate Performance
PDOP = get_PDOP_vec_sphere(X_ECI,T_vec,30,30,GMST0,15);

fit = (PDOP.'*PDOP)/length(PDOP) + max(PDOP)^2;
end