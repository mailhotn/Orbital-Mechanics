function [oeM] = osc2meNum(oeOsc,primary)
%osc2meNum returns the mean elements numerically by averaging the
%osculating over one orbital period
if nargin < 2
    primary = earth();
end
tN = 2*pi*sqrt(oeOsc(1)^3/primary.mu)*(1-1.5*primary.J2*(primary.Re/oeOsc(1))^2*(3-4*sind(oeOsc(3))^2));
tSpan = [0:10:tN];

Sat = SingleSat(oeOsc,primary);
Prop = Propagator(Sat);
[tVec,xOsc] = Prop.PropOeOsc3(tSpan);
M0 = xOsc(:,6) - 180/pi*tVec.*sqrt(primary.mu./xOsc(:,1).^3);
oeM = trapz(tVec,[xOsc(:,1:5),M0])/(tVec(end)-tVec(1)); % How to Average Mean anomaly???
oeM = oeM.';

