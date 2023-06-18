%% Define Satellite & Noise Model
oe = [7000, 0.01, 55, 10, 10, 10];
Sat = SingleSat(oe);
Prop = Propagator(Sat);

sigR = 5e-3;
sigV = 2e-5;
covEci = diag([sigR*ones(1,3),sigV*ones(1,3)]);
%% Prop & Contaminate
t = 0:10:86400;

eciTrue = PropEciJ2(t);
oeTrue = eci2oe(eciTrue);



