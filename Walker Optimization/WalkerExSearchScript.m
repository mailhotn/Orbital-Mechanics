%% Initialization
% Fitness Function of PDOP 'Integral'  'Mean'  'Max'
fitFunc = 'Mean';
maxPdop = 1000;
% Simulation Parameters
timeVector = 0:100:86164;
lonGs = 0;
elevMin  = 5;
relTol = 1e-6;
absTol = 1e-6;
% Genome Definition:         [P, F]


datafolder = 'C:\Users\User\Dropbox\Walker Optimization Data';

%% Massive For Loop
maxSats = 80;
minSats = 20;
for latGs = 10:10:80
    for T = minSats:maxSats
        GAsol = []; fit = []; GaRgtSol = [];
        % Optimize
        try
            % Run GA
            tic
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