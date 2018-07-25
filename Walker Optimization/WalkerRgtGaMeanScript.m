%% Initialization
% Fitness Function of PDOP 'Integral'  'Mean'  'Max'
fitFunc = 'Max';
maxPdop = 1000;
% Simulation Parameters
timeVector = 0:100:86164;
lonGs = 0;
elevMin  = 5;
relTol = 1e-6;
absTol = 1e-6;
% Genome Definition:         [P, F, j, inc, GMST0]
nParams = 5;
% Integer Values
intCon = [1,2,3]; % P,F,j
% GA Options
Options = gaoptimset('UseParallel',true,...
                     'PopulationSize',100);
%     'PlotFcns',{@gaplotscorediversity,@gaplotrange,@gaplotexpectation});
% Bounds
P_LB = 1;
% P_UB = floor(T/2); % set at each step
F_LB = 0;
% F_UB = floor(T/2)-1; % set at each step
inc_LB = 0;
inc_UB = 90;
j_LB = 14;
j_UB = 15;
GMST_LB = 0;
GMST_UB = 90;

datafolder = 'C:\Users\User\Dropbox\Walker Optimization Data';

%% Massive For Loop
maxSats = 80;
minSats = 20;
primeList = primes(maxSats);
for latGs = 10:10:80
    for T = minSats:maxSats
        if ~any(T==primeList)
            GAsol = []; fit = []; GaRgtSol = []; %#ok<NASGU>
            % Set T-dependant bounds
            P_UB = floor(T/2);
            F_UB = floor(T/2)-1;
            LB = [P_LB, F_LB, j_LB, inc_LB, GMST_LB];
            UB = [P_UB, F_UB, j_UB, inc_UB, GMST_UB];
            % Optimize
            try
                % Run GA
                tic
                [GAsol, fit, ~] = ga(@(x)WalkerFitnessRgtMean(x,T,timeVector,...
                    latGs,lonGs,elevMin,relTol,absTol,maxPdop,fitFunc),nParams,[],[],[],[],...
                    LB,UB,@(x)WalkerGaNonLinearConstraints(x,T),intCon,Options);
                optTime = toc;
                % Verbose Output Message
                c = clock;
                sma = CalcRgtElement([],0,GAsol(4),GAsol(3),1);
                disp([newline num2str(c(3)) '/' num2str(c(2)) ' ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6),2)...
                    newline 'Optimization Complete for T = ' num2str(T) ...
                    ' and Latitude = ' num2str(latGs)...
                    newline 'Fitness(' fitFunc '): ' num2str(fit)...
                    newline 'Elapsed Time: ' num2str((optTime-mod(optTime,60))/60) ' min '...
                    num2str(mod(optTime,60),3) ' sec'...
                    newline 'Optimal Solution i:T/P/F: '...
                    num2str(GAsol(4)) '°'  ':' num2str(T) '/' num2str(GAsol(1))...
                    '/' num2str(GAsol(2)) ...
                    newline 'Semimajor Axis: ' num2str(sma) ' km' ...
                    newline 'Earth Repeats: ' num2str(GAsol(3))])
                % Save Data
                GaRgtSol.nSatsT   = T;
                GaRgtSol.nPlanesP = GAsol(1);
                GaRgtSol.phasingF = GAsol(2);
                GaRgtSol.jRepeats = GAsol(3);
                GaRgtSol.inc      = GAsol(4);
                GaRgtSol.gmst0    = GAsol(5);
                GaRgtSol.fit      = fit;
                GaRgtSol.latGs    = latGs;
                GaRgtSol.lonGs    = lonGs;
                GaRgtSol.elevMin  = elevMin;
                save([datafolder '\WalkerMeanRgtSolLat_' num2str(latGs)...
                    '_T_' num2str(T) '.mat'],'GaRgtSol');
            catch ME
                % Ga Error!
                c = clock;
                fileID = fopen([datafolder '\error_'...
                    num2str(c(3)) '-' num2str(c(2)) '_' num2str(c(4)) ...
                    '.' num2str(c(5)) '.txt'],'w'); % create error file for IFTT phone notification
                fprintf(fileID,ME.message);
                rethrow(ME)
            end
        end
    end
end