function [ th ] = me2ta( Me , e, tol)
%me2ta Calculates the True anomaly (in degrees) corresponding 
% to the given Mean anomaly (in degrees)
if nargin < 3
    tol = 1e-14;
end
Me = wrapTo2Pi(Me*pi/180);

E = (Me < pi).*(Me+e) + (Me >= pi).*(Me-e);
dE = inf*ones(1,length(Me));

while norm(dE,inf) > tol
    dE = -(E-e.*sin(E)-Me)./(1-e.*cos(E));
    E  = E + dE;
end

th = 2*atan(sqrt((1+e)./(1-e)).*tan(E/2));
th = wrapTo360(th*180/pi);
end