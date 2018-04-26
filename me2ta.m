function [ th ] = me2ta( Me , e, tol)
%mean2true_anomaly Calculates the True anomaly (in radians) corresponding 
% to the given Mean anomaly (in radians)
if nargin < 3
    tol = 1e-14;
end
if Me < pi
    E = Me + e;
else
    E = Me - e;
end
dE = inf;
while abs(dE) > tol
    dE = -(E-e*sin(E)-Me)/(1-e*cos(E));
    E  = E + dE;
end
th = 2*atan(sqrt((1+e)/(1-e))*tan(E/2));
if th < 0
    th = 2*pi + th;
end
end