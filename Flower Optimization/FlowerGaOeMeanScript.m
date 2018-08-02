%% Initialization
% Fitness Function of PDOP 'Integral'  'Mean'  'Max'
PropParams.fitFunc = 'Integral';
PropParams.maxPdop = 1000;
% Simulation Parameters
PropParams.lonGs = 0;
PropParams.elevMin  = 5;
PropParams.relTol = 1e-6;
PropParams.absTol = 1e-6;
% Genome(9): nPetals-fN-fD-fH-w-i-altP-raan0-M0
nParams = 9;
% Integer Values
intCon = [1,2,3,4]; % nPetals-fN-fD-fH
% GA Options
Options = gaoptimset('UseParallel',true,...
                     'PopulationSize',100);
%     'PlotFcns',{@gaplotscorediversity,@gaplotrange,@gaplotexpectation});
% Bounds        nP,  fN,  fD,  fH,   w, inc, altP, raan0,  M0 
lowerBounds = [  1,   1,   1,   0,   0,   0,  500,     0,   0];
upperBounds = [nan, nan, nan, nan, 360,  90, 1000,   180, 180];

datafolder = 'C:\Users\User\Dropbox\Flower Optimization Data';

%% Massive For Loop
maxSats = 80;
minSats = 20;
primeList = primes(maxSats);
for latGs = 10:10:80
    for nSats = minSats:maxSats
        if ~any(nSats==primeList)
            GAsol = []; fit = []; GaRgtSol = []; %#ok<NASGU>
            % Set nSats-dependant bounds
            upperBounds(1:4) = [14*nSats, floor(nSats/2), floor(nSats/2), nSats-1];
            % Optimize
            try
                % Run GA
                tic
                [GAsol, fit, ~] = ga(@(x)WalkerFitnessRgtMean(x,nSats,PropParams),...
                    nParams,[],[],[],[],lowerBounds,upperBounds,...
                    @(x)FlowerGaNonLinearConstraints(x,nSats),intCon,Options);
                optTime = toc;
                % Verbose Output Message
                c = clock;
                sma = CalcRgtElement([],0,GAsol(4),GAsol(3),1);
                disp([newline num2str(c(3)) '/' num2str(c(2)) ' ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6),2)...
                    newline 'Optimization Complete for T = ' num2str(nSats) ...
                    ' and Latitude = ' num2str(latGs)...
                    newline 'Fitness(' fitFunc '): ' num2str(fit)...
                    newline 'Elapsed Time: ' num2str((optTime-mod(optTime,60))/60) ' min '...
                    num2str(mod(optTime,60),3) ' sec'...
                    newline 'Optimal Solution i:T/P/F: '...
                    num2str(GAsol(4)) '°'  ':' num2str(nSats) '/' num2str(GAsol(1))...
                    '/' num2str(GAsol(2)) ...
                    newline 'Semimajor Axis: ' num2str(sma) ' km' ...
                    newline 'Earth Repeats: ' num2str(GAsol(3))])
                % Save Data
                GaRgtSol.nSatsT   = nSats;
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
                    '_T_' num2str(nSats) '.mat'],'GaRgtSol');
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