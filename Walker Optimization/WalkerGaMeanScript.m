%% Initialization
% Simulation Parameters
Time = 0:600:5*86400;
latGs = 30;
lonGs = 0;
elevMin  = 15;
% Genome Definition:         [P, F, inc, alt, GMST0]
nParams = 5;
% Integer Values
intCon = [1,2]; % P,F
% GA Options
Options = gaoptimset('UseParallel',true...
                    ,'PopulationSize',50);
%     ,'PlotFcns',{@gaplotbestf,@gaplotscorediversity,@gaplotrange,@gaplotexpectation}...


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
delete([datafolder '\error.txt']) % clear IFTT error file
N = 300;
p = primes(N);
for T = 279:N
    if ~any(T == p)
        GAsol = []; fit = []; sol = []; %#ok<NASGU>
        % Set T-dependant bounds
        P_UB = floor(T/2);
        F_UB = floor(T/2)-1;
        LB = [P_LB, F_LB, inc_LB, alt_LB, GMST_LB];
        UB = [P_UB, F_UB, inc_UB, alt_UB, GMST_UB];
        % Optimize
        try
            tic
            [GAsol, fit, ~] = ga(@(x)WalkerFitnessMean(x,T,Time,latGs,lonGs,elevMin)...
                ,nParams,[],[],[],[],LB,UB,@(x)WalkerGaNonLinearConstraints(x,T),...
                intCon,Options);
            optTime = toc;
            c = clock;
            disp([newline num2str(c(3)) '/' num2str(c(2)) ' ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6),2)...
                newline 'Optimization Complete for T = ' num2str(T)...
                newline 'Maximum PDOP(fitness): ' num2str(fit)...
                newline 'Elapsed Time: ' num2str((optTime-mod(optTime,60))/60) ' min '...
                num2str(mod(optTime,60),3) ' sec'...
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
            save([datafolder '\Walker_Mean_Sol_T_' num2str(T) '.mat'],'sol');
        catch ME
            c = clock;
            fileID = fopen([datafolder '\error'...
                num2str(c(3)) '_' num2str(c(2)) '_' num2str(c(4)) '_' num2str(c(5)) '.txt'],'w'); % create error file for IFTT phone notification
            fprintf(fileID,ME.message);
            rethrow(ME)
        end
    end
end