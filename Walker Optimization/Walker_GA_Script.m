%% Optimization

T = 200;
% P = 20;
% Genome Definition:
% [P, F, inc, alt, GMST0]
nParams = 4;
incLB = 0;
incUB = 180;
altLB = 400;
altUB = 2000;

GMSTLB = 0;
GMSTUB = 60;
LB = [incLB, altLB, GMSTLB];
UB = [incUB, altUB, GMSTUB];

F = 2;

options = gaoptimset('UseParallel',true,'PlotFcns',{@gaplotbestf,@gaplotbestindiv}...
    ,'CrossoverFraction',0.8,'PopulationSize',100);
[GAsol, fit] = ga(@(x)WalkerFitness_sphere(x,T,P,F),nParams,[],[],[],[],LB,UB,[],[],options);