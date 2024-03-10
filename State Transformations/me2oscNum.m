function [oeOsc1,fVal] = me2oscNum(oeM)
%me2oscNum finds the oscullating elements at time 0, that result in the
% desired mean elements, based on a simulation from time -T to time 0

% Convoluted way of saying this is a numerical mean to Oscullating
% transformation
primary = earth();
boundVec = oeM*primary.J2;
boundVec(2) = primary.J2;
boundVec(3:end) = 1;
ub = oeM +  3*boundVec;
lb = oeM - 3*boundVec;
lb(2:3) = 0;
options = optimoptions('fmincon',Display='off',...
    OptimalityTolerance=1e-8,...
    Algorithm='interior-point');
[oeOsc0,fVal] = fmincon(@(x)meanSim(x,oeM),oeM,[],[],[],[],lb,ub,[],options);
oeOsc1 = oeOsc0;
% add half orbit ?
Sat = SingleSat(oeOsc0,primary);
Prop = Propagator(Sat);
primary = earth();
tN = 2*pi*sqrt(oeOsc0(1)^3/primary.mu)*(1-1.5*primary.J2*(primary.Re/oeOsc0(1))^2*(3-4*sind(oeOsc0(3))^2));
tVec = [0:10:tN];
[~,xOsc] = Prop.PropOeOsc3(tVec);
M0 = xOsc(end,6) - 180/pi*tVec(end).*sqrt(primary.mu./xOsc(end,1).^3);
oeOsc1 = [xOsc(end,1:5),M0].';


% oeOsc1 = oeOsc0 + one half Orbit

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