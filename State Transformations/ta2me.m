function [ Me ] = ta2me(X1, X2)
%true2mean_anomaly calculates the Mean anomaly of a spacecraft given it's
%true anomaly and eccentricity.  All I/O angles are degrees.
if nargin < 2
    X2 = [];
end

if ~isempty(X2) % X1 is th, X2 is e
    th = X1;
    e = X2;
else % no X2 => X1 is full state
    if size(X1,1) ~= 6 && size(X1,2) == 6
        X1 = X1.';
    end
    th = X1(6,:);
    e = X1(2,:);
end

th = wrapTo2Pi(th*pi/180);

E = 2*atan(sqrt((1-e)./(1+e)).*tan(th/2));
M = E - e.*sin(E);

M = wrapTo360(M*180/pi);
if ~isempty(X2)
    Me = M;
else
    Me = [X1(1:5,:); M];
end