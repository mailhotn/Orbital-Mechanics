function [ th ] = me2ta( X1 , X2, tol)
%me2ta Calculates the True anomaly (in degrees) corresponding 
% to the given Mean anomaly (in degrees)
maxIter = 20;
if nargin < 2
    X2 = [];
end
if nargin < 3
    tol = 1e-14;
end
if ~isempty(X2) % X1 is M, X2 is e
    M = X1;
    e = X2;
else % no X2 => X1 is full state
    if size(X1,1) ~= 6 && size(X1,2) == 6
        X1 = X1.';
    end
    M = X1(6,:);
    e = X1(2,:);
end

M = wrapTo2Pi(M*pi/180);

E = (M < pi).*(M+e) + (M >= pi).*(M-e);
dE = inf(1,length(M));
iter = 0;
while (norm(dE,inf) > tol) && (iter < maxIter)
    dE = -(E-e.*sin(E)-M)./(1-e.*cos(E));
    E  = E + dE;
    iter = iter + 1;
end

E = wrapTo2Pi(E);
th = 2*atan(sqrt((1+e)./(1-e)).*tan(E/2));
th = wrapTo360(th*180/pi);
if ~isempty(X2)
    th = th;
else
    th = [X1(1:5,:); th];
end