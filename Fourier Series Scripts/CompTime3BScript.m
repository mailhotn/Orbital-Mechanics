clear
%% Setup Molniya
primary = Earth;
third = {Moon,Sun};
% T = 86164/2;
% sma = ((T/2/pi)^2*primary.mu)^(1/3);
% ic = [sma, 0.74, 63.4, 30, 270, 10].'; % Molniya

T = 86164;
sma = ((T/2/pi)^2*primary.mu)^(1/3);
ic = [sma, 0.001, 5, 30, 30, 10].'; % GEO

Sat = SingleSat(ic,primary,third);
Prop = Propagator(Sat);
t = 0:1800:86400*30;
kMax = 4;

nTest = 100;
numTime = 0;
fourTime = 0;
singleTime = 0;
%% Propagate Molniya
for iTest = 1:nTest
% Prop using Good solution
tic
[~,X1] = Prop.PropEci3B(t,'Stable');
oeN = eci2oe(X1.',[],primary,'me');
oeN(6,:) = unwrap(oeN(6,:)*pi/180)*180/pi - sqrt(primary.mu./oeN(1,:).^3).*t*180/pi;
pTime = toc;
numTime = numTime + pTime;

% Prop using Fourier solution
tic
[~,oeF] = Prop.PropOeFourier3B(t,kMax);
oeF = oeF.';
pTime = toc;
fourTime = fourTime + pTime;
% Prop using singly averaged
tic
[~,oeS] = Prop.PropOeMeanSingle3B(t);
oeS = oeS.';
pTime = toc;
singleTime = singleTime + pTime;
end
fourTime/numTime
singleTime/numTime 
numTime/fourTime
numTime/singleTime