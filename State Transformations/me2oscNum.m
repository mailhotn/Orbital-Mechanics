function [oeOsc1,fVal] = me2oscNum(oeM)
%me2oscNum finds the oscullating elements at time 0, that result in the
% desired mean elements, based on a simulation from time -T to time 0

% Convoluted way of saying this is a numerical mean to Oscullating
% transformation
primary = earth();
boundVec = oeM*primary.J2;
boundVec(2) = primary.J2;
ub = oeM +  boundVec;
lb = oeM - boundVec;
lb(2:3) = 0;
[oeOsc0,fVal] = fmincon(@(x)meanSim(x,oeM),oeM,[],[],[],[],lb,ub);
oeOsc1 = oeOsc0;
% Sat = SingleSat(oeOsc0,primary);
% Prop = Propagator(Sat);


% oeOsc1 = oeOsc0 + one Orbit

end

function [oeError] = meanSim(x,oeM)
% simulate 1 orbit with proposed elements, average elements, output square
% error between average and target mean

primary = earth();
tN = 2*pi*sqrt(x(1)^3/primary.mu)*(1-1.5*primary.J2*(primary.Re/x(1))^2*(3-4*sind(x(3))^2));
tSpan = [0:10:tN];

Sat = SingleSat(x,primary);
Prop = Propagator(Sat);
[tVec,xOsc] = Prop.PropOeOsc3(tSpan);
M0 = xOsc(:,6) - 180/pi*tVec.*sqrt(primary.mu./xOsc(:,1).^3);
xM = trapz(tVec,[xOsc(:,1:5),M0])/(tVec(end)-tVec(1)); % How to Average Mean anomaly???

unitVec = [1/primary.Re,1,ones(1,4)*pi/180]; 
oeError = ((xM-oeM.').*unitVec)*((xM.'-oeM).*unitVec.')*1000; 

end