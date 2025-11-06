function [oeM] = osc2meNum(oeOsc,primary)
%osc2meNum returns the mean elements numerically by averaging the
%osculating over one orbital period
if nargin < 2
    primary = earth();
end
% tN = 2*pi*sqrt(oeOsc(1)^3/primary.mu)*(1-1.5*primary.J2*(primary.Re/oeOsc(1))^2*(3-4*sind(oeOsc(3))^2));
T = 2*pi*sqrt(oeOsc(1)^3/primary.mu);
% nT = floor(T/100);
tSpan = [0,T];

Sat = SingleSat(oeOsc,primary);
Prop = Propagator(Sat,[],1e-5,1e-6);
[tVec,xOsc] = Prop.PropOeOsc3(tSpan);
xM = nan(1,6);
% Unwrap O, w, M before averaging
xOsc(:,4:6) = 180/pi*unwrap(pi/180*xOsc(:,4:6));
% Use mean n to get M0
xM(1:5) = trapz(tVec,xOsc(:,1:5))/(tVec(end)-tVec(1));
M0 = xOsc(:,6) - 180/pi*tVec.*sqrt(primary.mu./xM(1).^3); %M0 = Mosc - nMt - mean mean motion
xM(6) = trapz(tVec,M0)/(tVec(end)-tVec(1));
% Don't use mean n to get M0
% M0 = xOsc(:,6) - 180/pi*tVec.*sqrt(primary.mu./xOsc(:,1).^3); %M0 = Mosc - nt
% xM = trapz(tVec,[xOsc(:,1:5),M0])/(tVec(end)-tVec(1)); % How to Average Mean anomaly???
xM(4:5) = wrapTo360(xM(4:5));
oeM = xM.';

