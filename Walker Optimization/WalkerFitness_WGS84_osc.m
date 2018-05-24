function [ fit ] = WalkerFitness_WGS84_osc(x, T, Time, lat_gs, lon_gs, e_min, reltol, abstol)
%WalkerFitness_sphere Simulates a Walker Constellation and calculates the
%fitness

% Initialization
P = x(1);
F = x(2);
inclination = x(3);
altitude = x(4);
GMST0 = x(5);

WC = WalkerConstellation(T,P,F,inclination,altitude);
Prop = Propagator(WC,reltol,abstol);

% Propagate
[T_vec,X_ECI] = Prop.prop_ECI_J2(Time);

% Evaluate Performance
PDOP = get_PDOP_vec_WGS84(X_ECI,T_vec,lat_gs,lon_gs,GMST0,e_min);
fit = max(PDOP);
end