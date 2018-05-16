%% Optimization

T = 200;
% Genome Definition:         [P, F, inc, alt, GMST0]
% 
nParams = 5;
% Bounds
P_LB = 1;
P_UB = floor(T/2);
F_LB = 0;
F_UB = floor(T/2)-1;
inc_LB = 0;
inc_UB = 180;
alt_LB = 400;
alt_UB = 2000;
GMST_LB = 0;
GMST_UB = 60;
LB = [P_LB, F_LB, inc_LB, alt_LB, GMST_LB];
UB = [P_UB, F_UB, inc_UB, alt_UB, GMST_UB];
% Integer Values
intcon = [1,2];
% Options
options = gaoptimset('UseParallel',true,'PlotFcns',{@gaplotbestf,@gaplotbestindiv}...
    ,'CrossoverFraction',0.8,'PopulationSize',100);
% Optimize
[GAsol, fit] = ga(@(x)WalkerFitness_sphere(x,T),nParams,[],[],[],[],...
                 LB,UB,@(x)Walker_GA_nonlcon(x,T),intcon,options);