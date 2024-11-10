function [oeOsc1,fVal,exitflag,output] = me2oscNum(oeM)
%me2oscNum finds the oscullating elements at time 0, that result in the
% desired mean elements, based on a simulation from time -T to time 0

% Convoluted way of saying this is a numerical mean to Oscullating
% transformation for one satellite

primary = earth();
eta = sqrt(1-oeM(2)^2);
% set bounds
boundVec = 2*oeM*primary.J2;
boundVec(2) = 3*primary.J2;
boundVec(3:end) = [1,1,4,3];


A = [];
b = [];
% Alternate bounds for low ecc
if oeM(2) < 1e-2
    delLon = 4;
    boundVec(3:6) = [0.5,0.5,25,25];
    A = [0 0 0 0 1 1;
        0 0 0 0 -1 -1];
    b = [(oeM(5)+oeM(6)) + delLon;
        -(oeM(5)+oeM(6)) + delLon];
end
% Finish Bounds
ub = oeM + 3*boundVec/eta^4;
lb = oeM - 3*boundVec/eta^4;
lb(2:3) = lb(2:3).*(lb(2:3)>0); % bring lower bounds up from 0

options = optimoptions('fmincon',Display='off',...
    OptimalityTolerance=1e-8,...
    Algorithm='interior-point');
[oeOsc0,fVal,exitflag,output] = fmincon(@(x)meanSim(x,oeM),oeM,A,b,[],[],lb,ub,[],options);
oeOsc0(4:5) = wrapTo360(oeOsc0(4:5));
oeOsc1 = oeOsc0;
% add half orbit ?
% Sat = SingleSat(oeOsc0,primary);
% Prop = Propagator(Sat);
% primary = earth();
% % tN = 2*pi*sqrt(oeOsc0(1)^3/primary.mu)*(1-1.5*primary.J2*(primary.Re/oeOsc0(1))^2*(3-4*sind(oeOsc0(3))^2)); % nodal period
% T = 2*pi*sqrt(oeOsc0(1)^3/primary.mu);
% % tVec = [0:10:tN];
% nT = floor(T/10);
% tVec = linspace(0,T,nT);
% [~,xOsc] = Prop.PropOeOsc3(tVec);
% M0 = xOsc(end,6) - 180/pi*tVec(end).*sqrt(primary.mu./xOsc(end,1).^3);
% oeOsc1 = [xOsc(end,1:5),M0].';


% oeOsc1 = oeOsc0 + one Orbit

end

function [oeError] = meanSim(x,oeM)
% simulate 1 orbit with proposed elements, average elements, output square
% error between average and target mean

primary = earth();
% tN = 2*pi*sqrt(x(1)^3/primary.mu)*(1-1.5*primary.J2*(primary.Re/x(1))^2*(3-4*sind(x(3))^2));
T = 2*pi*sqrt(x(1)^3/primary.mu);
nT = floor(T/100);
tSpan = linspace(0,T,nT);

Sat = SingleSat(x,primary);
Prop = Propagator(Sat);
[tVec,xOsc] = Prop.PropOeOsc3(tSpan);
xM = nan(1,6);

% Fix edge cases resulting from wrapping, where x>360,x<0; x = AOP, RAAN
if x(4) > 360 && xOsc(1,4) < 360
    xOsc(:,4) = xOsc(:,4) + 360;
elseif x(4) < 0 && xOsc(1,4) > 0
    xOsc(:,4) = xOsc(:,4) - 360;
end
if x(5) > 360 && xOsc(1,5) < 360
    xOsc(:,5) = xOsc(:,5) + 360;
elseif x(5) < 0 && xOsc(1,5) > 0
    xOsc(:,5) = xOsc(:,5) - 360;
end

% Unwrap O, w, M before averaging
xOsc(:,4:6) = 180/pi*unwrap(pi/180*xOsc(:,4:6));

% Use mean n to get M0
xM(1:5) = trapz(tVec,xOsc(:,1:5))/(tVec(end)-tVec(1));
M0 = xOsc(:,6) - 180/pi*tVec.*sqrt(primary.mu./xM(1).^3); %M0 = Mosc - nMt : mean mean motion
xM(6) = trapz(tVec,M0)/(tVec(end)-tVec(1));

% Don't use mean n to get M0
% M0 = xOsc(:,6) - 180/pi*tVec.*sqrt(primary.mu./xOsc(:,1).^3); %M0 = Mosc - nt
% xM = trapz(tVec,[xOsc(:,1:5),M0])/(tVec(end)-tVec(1)); % How to Average Mean anomaly???

% Make sure AOP & RAAN are continuous at 0,360 ????
% if xM(4) > 359
%     xM(4) = wrapTo180(xM(4));
% end
% if xM(5) > 360 % Seemingly important for edgecase AOP = 0
%     xM(5) = wrapTo180(xM(5));
% end
unitVec = [1/primary.Re,1,ones(1,4)*pi/180]; 
oeError = ((xM-oeM.').*unitVec)*((xM.'-oeM).*unitVec.')*1000; 

end