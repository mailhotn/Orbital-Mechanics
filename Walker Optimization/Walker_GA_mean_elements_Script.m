%% Initialization
% Simulation Parameters
Time = 0:600:5*86400;
lat_gs = 30;
lon_gs = 0;
e_min  = 15;
% Genome Definition:         [P, F, inc, alt, GMST0]
nParams = 5;
% Integer Values
intCon = [1,2]; % P,F
% GA Options
options = gaoptimset('UseParallel',true);%...
%     ,'PlotFcns',{@gaplotbestf,@gaplotscorediversity,@gaplotrange,@gaplotexpectation}...
%     ,'PopulationSize',1000);

% Bounds
P_LB = 1;
% P_UB = floor(T/2); % set at each step
F_LB = 0;
% F_UB = floor(T/2)-1; % set at each step
inc_LB = 0;
inc_UB = 170;
alt_LB = 400;
alt_UB = 2000;
GMST_LB = 0;
GMST_UB = 60;

datafolder = 'C:\Users\User\Dropbox\Walker Optimization Data';

%% Massive For Loop
N = 300;
p = primes(N);
for T = 40:N
    if ~any(T == p)
        GAsol = []; fit = []; sol = []; %#ok<NASGU>
        % Set T-dependant bounds
        P_UB = floor(T/2);
        F_UB = floor(T/2)-1;
        LB = [P_LB, F_LB, inc_LB, alt_LB, GMST_LB];
        UB = [P_UB, F_UB, inc_UB, alt_UB, GMST_UB];
        % Optimize
        tic
        [GAsol, fit, ~] = ga(@(x)WalkerFitness_WGS84_mean(x,T,Time,lat_gs,lon_gs,e_min)...
            ,nParams,[],[],[],[],LB,UB,@(x)Walker_GA_nonlcon(x,T),intCon,options);
        opt_time = toc;
        disp([newline 'Optimization Complete for T = ' num2str(T)...
            newline 'Elapsed Time: ' num2str(opt_time) ' seconds'...
            newline 'Optimal Solution (i/alt:T/P/F): '...
            num2str(GAsol(3)) '°/' num2str(GAsol(4)) ':' num2str(T) '/' num2str(GAsol(1))...
            '/' num2str(GAsol(2)) newline])
        sol.T = T;
        sol.P = GAsol(1);
        sol.F = GAsol(2);
        sol.inc = GAsol(3);
        sol.alt = GAsol(4);
        sol.GMST0 = GAsol(5);
        sol.fit = fit;
        save([datafolder '\Walker_Mean_Sol_T_' num2str(T) '_fit_' num2str(fit) '.mat'],'sol');
    end
end