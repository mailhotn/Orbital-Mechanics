function [ th ] = me2ta( Me , e, tol)
%me2ta Calculates the True anomaly (in degrees) corresponding 
% to the given Mean anomaly (in degrees)
maxIter = 20;
if nargin < 3
    tol = 1e-14;
end
Me = wrapTo2Pi(Me*pi/180);

E = (Me < pi).*(Me+e) + (Me >= pi).*(Me-e);
dE = inf(1,length(Me));
iter = 0;
while (norm(dE,inf) > tol) && (iter < maxIter)
    dE = -(E-e.*sin(E)-Me)./(1-e.*cos(E));
    E  = E + dE;
    iter = iter + 1;
end

E = wrapTo2Pi(E);
th = 2*atan(sqrt((1+e)./(1-e)).*tan(E/2));
th = wrapTo360(th*180/pi);
end