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
primeList = primes(maxSats);
T         = zeros(1,maxSats-minSats+1);
P         = zeros(1,maxSats-minSats+1);
F         = zeros(1,maxSats-minSats+1);
inc       = zeros(1,maxSats-minSats+1);
alt       = zeros(1,maxSats-minSats+1);
GMST0     = zeros(1,maxSats-minSats+1);
meanFit   = zeros(1,maxSats-minSats+1);
oscFit    = zeros(1,maxSats-minSats+1);
oscFitPat = zeros(1,maxSats-minSats+1);

% Pattern Search Options
incLB = 0;
incUB = 170;
altLB = 400;
altUB = 2000;
gmstLB = 0;
gmstUB = 60;
LB = [incLB, altLB, gmstLB];
UB = [incUB, altUB, gmstUB];
Options = psoptimset('UseParallel',true,'PlotFcn',{@psplotbestf});

datafolder = 'C:\Users\User\Dropbox\Walker Optimization Data';

%% Go over all solutions & sim Osculating to get new fitness
% Then use patternsearch with osc elements to improve fitness.
for iSats = (minSats+1):(maxSats+1)
    if ~any((iSats-1) == primeList)
        load([datafolder '\Walker_Mean_Sol_T_' num2str(iSats-1) '.mat']);
        if sol.fit == 20
            T(iSats-minSats)         = sol.T;
            P(iSats-minSats)         = nan;
            F(iSats-minSats)         = nan;
            inc(iSats-minSats)       = nan;
            alt(iSats-minSats)       = nan;
            GMST0(iSats-minSats)     = nan;
            meanFit(iSats-minSats)   = nan;
            oscFit(iSats-minSats)    = nan;
            oscFitPat(iSats-minSats) = nan;
        else
            T(iSats-minSats)        = sol.T;
            P(iSats-minSats)        = sol.P;
            F(iSats-minSats)        = sol.F;
            inc(iSats-minSats)      = sol.inc;
            alt(iSats-minSats)      = sol.alt;
            GMST0(iSats-minSats)    = sol.GMST0;
            meanFit(iSats-minSats)  = sol.fit;
            oscFit(iSats-minSats)   = WalkerFitness_WGS84_osc([P(iSats-minSats),F(iSats-minSats),inc(iSats-minSats)...
                ,alt(iSats-minSats),GMST0(iSats-minSats)],T(iSats-minSats),timeVector,latGs,lonGs,elevMin,relTol,absTol);
            initialGuessPat = [sol.inc, sol.alt, sol.GMST0];
            if ~exist([datafolder '\GA Pattern Data\WalkerOscPatSol_T_'...
                    num2str(sol.T) '.mat'],'file')
                try
                    tic
                    [xPat, fitPat]  = patternsearch(@(x)WalkerFitnessOscPat(x, sol.T, sol.P,...
                        sol.F, timeVector, latGs, lonGs, elevMin, relTol, absTol),...
                        initialGuessPat,[],[],[],[],LB,UB,[],Options);
                    optTime = toc;
                    oscFitPat(iSats-minSats) = fitPat;
                    solPat.T     = sol.T;
                    solPat.P     = sol.P;
                    solPat.F     = sol.F;
                    solPat.inc   = xPat(1);
                    solPat.alt   = xPat(2);
                    solPat.GMST0 = xPat(3);
                    solPat.fit   = fitPat;
                    save([datafolder '\GA Pattern Data\WalkerOscPatSol_T_' num2str(solPat.T) '.mat'],'solPat');
                    c = clock;
                    disp([newline num2str(c(3)) '/' num2str(c(2)) ' ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6),2)...
                        newline 'Pattern Search Complete for T = ' num2str(sol.T)...
                        newline 'Mean Elements Maximum PDOP: ' num2str(sol.fit)...
                        newline 'Osculating Elements Maximum PDOP: ' num2str(oscFit(iSats-minSats))...
                        newline 'Pattern Search Maximum PDOP: ' num2str(fitPat)...
                        newline 'Elapsed Time: ' num2str((optTime-mod(optTime,60))/60) ' min '...
                        num2str(mod(optTime,60),3) ' sec'...
                        newline 'Optimal Solution (i/alt:T/P/F): '...
                        num2str(solPat.inc) '°/' num2str(solPat.alt) ':' num2str(solPat.T)...
                        '/' num2str(solPat.P) '/' num2str(solPat.F) newline])
                catch ME
                    c = clock;
                    fileID = fopen([datafolder '\GA Pattern Data\error'...
                        num2str(c(3)) '_' num2str(c(2)) '_' num2str(c(4)) '_' num2str(c(5)) '.txt'],'w'); % create error file for IFTT phone notification
                    fprintf(fileID,ME.message);
                    rethrow(ME)
                end
            else
                 load([datafolder '\GA Pattern Data\WalkerOscPatSol_T_'...
                     num2str(sol.T) '.mat']);
                 oscFitPat(iSats-minSats) = solPat.fit;
                 disp(['Pattern Search Already Completed for T = ' num2str(sol.T)...
                        newline 'Mean Elements Maximum PDOP: ' num2str(sol.fit)...
                        newline 'Osculating Elements Maximum PDOP: ' num2str(oscFit(iSats-minSats))...
                        newline 'Pattern Search Maximum PDOP: ' num2str(solPat.fit)...
                        newline 'Optimal Solution (i/alt:T/P/F): '...
                        num2str(solPat.inc) '°/' num2str(solPat.alt) ':' num2str(solPat.T)...
                        '/' num2str(solPat.P) '/' num2str(solPat.F) newline])
            end
        end
    else
        T(iSats-minSats)        = iSats-1;
        P(iSats-minSats)        = nan;
        F(iSats-minSats)        = nan;
        inc(iSats-minSats)      = nan;
        alt(iSats-minSats)      = nan;
        GMST0(iSats-minSats)    = nan;
        meanFit(iSats-minSats)  = nan;
        oscFit(iSats-minSats)   = nan;
    end
end
indBad = find(oscFit == 20);
P(indBad)        = nan;
F(indBad)        = nan;
inc(indBad)      = nan;
alt(indBad)      = nan;
GMST0(indBad)    = nan;
meanFit(indBad)  = nan;
oscFit(indBad)   = nan;
%% Plot Results
%Fitness Plot
figure()
semilogy(T,meanFit,'o',T,oscFit,'o',T,oscFitPat,'s');
xlabel('# Satellites')
ylabel('Maximum PDOP')
legend('Mean Sim','Osculating Sim','Hybrid Optimization')
% inclination
figure()
plot(T,inc,'o');
xlabel('# Satellites')
ylabel('Inclination [°]')
% Altitude
figure()
plot(T,alt,'o');
xlabel('# Satellites')
ylabel('Altitude [km]')