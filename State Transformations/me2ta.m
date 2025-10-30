function [ th ] = me2ta( X1 , X2, tol)
%me2ta Calculates the True anomaly (in degrees) corresponding 
% to the given Mean anomaly (in degrees)

% X1 = M
% X2 = e
% if no X2 is provided, X1 is instead the full state

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

M = mod(M*pi/180,2*pi);

E = (M < pi).*(M+e) + (M >= pi).*(M-e);
dE = inf(1,length(M));
iter = 0;
while (norm(dE,inf) > tol) && (iter < maxIter)
    dE = -(E-e.*sin(E)-M)./(1-e.*cos(E));
    E  = E + dE;
    iter = iter + 1;
end

E = mod(E,2*pi);
th = 2*atan(sqrt((1+e)./(1-e)).*tan(E/2));
th = mod(th*180/pi,360);
if ~isempty(X2)
    th = th;
else
    th = [X1(1:5,:); th];
end