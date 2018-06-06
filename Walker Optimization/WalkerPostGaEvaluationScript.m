%% Initializaton
% Simulation Parameters
timeVector = 0:10:5*86400;
latGs = 30;
lonGs = 0;
elevMin  = 15;
relTol = 1e-6;
absTol = 1e-6;

% Initialize data vectors
minSats = 50;
maxSats  = 300;
nCons = maxSats - minSats + 1;
primeList = primes(maxSats);
T         = nan(1,nCons);
P         = nan(1,nCons);
F         = nan(1,nCons);
gaInc     = nan(1,nCons);
gaAlt     = nan(1,nCons);
gaGmst0   = nan(1,nCons);
patInc    = nan(1,nCons);
patAlt    = nan(1,nCons);
patGmst0  = nan(1,nCons);
meanFit   = nan(1,nCons);
oscFit    = nan(1,nCons);
oscFitPat = nan(1,nCons);

% Pattern Search Options
incLB = 0;
incUB = 170;
altLB = 400;
altUB = 2000;
gmstLB = 0;
gmstUB = 60;
LB = [incLB, altLB, gmstLB];
UB = [incUB, altUB, gmstUB];
Options = psoptimset('UseParallel',false,'PlotFcn',{@psplotbestf});

datafolder = 'C:\Users\User\Dropbox\Walker Optimization Data';

%% Go over all solutions & sim Osculating to get new fitness
% Load previously run pattersearch solutions
for iSats = 1:nCons
    if ~any((iSats + minSats - 1) == primeList)
        load([datafolder '\Walker_Mean_Sol_T_' num2str(iSats + minSats -1) '.mat']);
        if sol.fit ~= 20
            % assign values to vectors
            T(iSats)        = sol.T;
            P(iSats)        = sol.P;
            F(iSats)        = sol.F;
            gaInc(iSats)    = sol.inc;
            gaAlt(iSats)    = sol.alt;
            gaGmst0(iSats)  = sol.GMST0;
            meanFit(iSats)  = sol.fit;
            % Osc Sim
            oscFit(iSats)   = WalkerFitnessOsc([P(iSats),...
                F(iSats),gaInc(iSats),gaAlt(iSats),...
                gaGmst0(iSats)],T(iSats),timeVector,...
                latGs,lonGs,elevMin,relTol,absTol);
            % check for pattern solution
            if exist([datafolder '\GA Pattern Data\WalkerOscPatSol_T_'...
                    num2str(T(iSats)) '.mat'],'file')
                % load pattern solution
                load([datafolder '\GA Pattern Data\WalkerOscPatSol_T_'...
                    num2str(T(iSats)) '.mat']);
                % assign values to vectors
                oscFitPat(iSats) = solPat.fit;
                patInc(iSats)    = solPat.inc;
                patAlt(iSats)    = solPat.alt;
                patGmst0(iSats)  = solPat.GMST0;
                % display output
                disp(['Pattern Search Already Completed for T = ' num2str(sol.T)...
                    newline 'Mean Elements Maximum PDOP: ' num2str(sol.fit)...
                    newline 'Osculating Elements Maximum PDOP: ' num2str(oscFit(iSats))...
                    newline 'Pattern Search Maximum PDOP: ' num2str(solPat.fit)...
                    newline 'Optimal Solution (i/alt:T/P/F): '...
                    num2str(solPat.inc) '°/' num2str(solPat.alt) ':' num2str(solPat.T)...
                    '/' num2str(solPat.P) '/' num2str(solPat.F) newline])
            end % end case of pattern already computed
        else % fitness = 20, ga failure
            T(iSats)         = sol.T;
        end % end if fit ~= 20
    else % prime, no data
        T(iSats)        = iSats-1;
    end  % end prime check
end % end for loop

%% Run Pattern Search for all cases where no solution exists yet
parfor iSats = 1:nCons
    patternSolExists = exist([datafolder '\GA Pattern Data\WalkerOscPatSol_T_'...
            num2str(T(iSats)) '.mat'],'file');
    isPrime = any((iSats + minSats - 1) == primeList);
    isBad = isnan(oscFit(iSats));
    if ~isPrime && ~patternSolExists && ~isBad
        try
            initialGuessPat = [gaInc(iSats), gaAlt(iSats), gaGmst0(iSats)];
            tic
            [xPat, fitPat]  = patternsearch(@(x)WalkerFitnessOscPat(x,...
                T(iSats), P(iSats), F(iSats), timeVector, latGs, lonGs,...
                elevMin, relTol, absTol), initialGuessPat, [], [], [], [],...
                LB, UB, [], Options);
            optTime = toc;
            oscFitPat(iSats) = fitPat;
            patInc(iSats)    = xPat(1);
            patAlt(iSats)    = xPat(2);
            patGmst0(iSats)  = xPat(3);
            c = clock;
            disp([newline num2str(c(3)) '/' num2str(c(2)) ' ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6),2)...
                newline 'Pattern Search Complete for T = ' num2str(T(iSats))...
                newline 'GA w/ Mean Elements Maximum PDOP: ' num2str(meanFit(iSats))...
                newline 'Osculating Elements Maximum PDOP: ' num2str(oscFit(iSats))...
                newline 'Pattern & GA Hybrid Maximum PDOP: ' num2str(oscFitPat(iSats))...
                newline 'Elapsed Time: ' num2str((optTime-mod(optTime,60))/60) ' min '...
                num2str(mod(optTime,60),3) ' sec'...
                newline 'Optimal Solution (i/alt:T/P/F): '...
                num2str(patInc(iSats)) '°/' num2str(patAlt(iSats))...
                ':' num2str(T(iSats)) '/' num2str(P(iSats)) '/'...
                num2str(F(iSats)) newline])
        catch ME
            c = clock;
            fileID = fopen([datafolder '\GA Pattern Data\error'...
                num2str(c(3)) '_' num2str(c(2)) '_' num2str(c(4)) '_' num2str(c(5)) '.txt'],'w'); % create error file for IFTT phone notification
            fprintf(fileID,ME.message);
            rethrow(ME)
        end % end try/catch
    end % end pattern search
end % end parfor

% Save Pattern Solutions
for iSats = 1:nCons
    patternSolExists = exist([datafolder '\GA Pattern Data\WalkerOscPatSol_T_'...
            num2str(T(iSats)) '.mat'],'file');
    isPrime = any((iSats + minSats - 1) == primeList);
    isBad = isnan(oscFit(iSats));
    if ~isPrime && ~patternSolExists && ~isBad
        solPat.fit   = oscFitPat(iSats);
        solPat.T     = T(iSats);
        solPat.P     = P(iSats);
        solPat.F     = F(iSats);
        solPat.inc   = patInc(iSats);
        solPat.alt   = patAlt(iSats);
        solPat.GMST0 = patGmst0(iSats);
        save([datafolder '\GA Pattern Data\WalkerOscPatSol_T_' num2str(solPat.T) '.mat'],'solPat');
    end
end
    
% indBad = find(oscFit == 20);
% P(indBad)        = nan;
% F(indBad)        = nan;
% gaInc(indBad)      = nan;
% gaAlt(indBad)      = nan;
% gaGmst0(indBad)    = nan;
% meanFit(indBad)  = nan;
% oscFit(indBad)   = nan;
%% Plot Results
%Fitness Plot
figure()
semilogy(T,meanFit,'o',...
         T,oscFit,'o',...
         T,oscFitPat,'s');
xlabel('# Satellites')
ylabel('Maximum PDOP')
legend('Mean Sim','Osculating Sim','Hybrid Optimization')
% inclination
figure()
plot(T,gaInc,'o',...
     T,patInc,'o');
xlabel('# Satellites')
ylabel('Inclination [°]')
legend('Mean Sim','Hybrid Optimization')
% Altitude
figure()
plot(T,gaAlt,'o',...
     T,patAlt,'o');
xlabel('# Satellites')
ylabel('Altitude [km]')
legend('Mean Sim','Hybrid Optimization')