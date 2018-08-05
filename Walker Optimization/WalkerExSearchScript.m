%% Initialization
PropParams.maxPdop = 1000;
PropParams.timeVec = 0:10:86164;
lonGs = 0;
PropParams.elevMin  = 5;
PropParams.relTol = 1e-6;
PropParams.absTol = 1e-6;
PropParams.datafolder = 'C:\Users\User\Dropbox\Walker Optimization Data';

latList = 10:10:80;
Rgt.jRepeats = 14;
Rgt.kDays    = 1;
nOrbits = 6;

%% Massive For Loop
maxSats = 100;
minSats = 20;
parfor iLat = 1:length(latList)
    GroundS = struct();
    GroundS.lat = latList(iLat);
    GroundS.lon = 0;
    ecc = 0;
    [inc, sma, raan0] = VirtualSatOe(ecc,Rgt,GroundS,nOrbits);
    for nSats = minSats:maxSats
        % Optimize
        try
            % Run Search for T & latGs
            tic
            ExSol = WalkerExSearch(nSats, inc, sma, raan0, latList(iLat), PropParams);
            optTime = toc;
            % Verbose Output Message
            c = clock;
            disp([newline num2str(c(3)) '/' num2str(c(2)) ' ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6),2)...
                newline 'Optimization Complete for T = ' num2str(nSats) ...
                ' and Latitude = ' num2str(latList(iLat))...
                newline 'Fitness: ' num2str(ExSol.fit)...
                newline 'Elapsed Time: ' num2str((optTime-mod(optTime,60))/60) ' min '...
                num2str(mod(optTime,60),3) ' sec'...
                newline 'Optimal Solution i:T/P/F: '...
                num2str(inc) '°'  ':' num2str(nSats) '/' num2str(ExSol.optP)...
                '/' num2str(ExSol.optF) ...
                newline 'Semimajor Axis: ' num2str(sma) ' km' ...
                newline 'Earth Repeats: ' num2str(Rgt.jRepeats)])
        catch ME
            % Error!
            c = clock;
            fileID = fopen([PropParams.datafolder '\error_'...
                num2str(c(3)) '-' num2str(c(2)) '_' num2str(c(4)) ...
                '.' num2str(c(5)) '.txt'],'w'); % create error file for IFTT phone notification
            fprintf(fileID,ME.message);
            rethrow(ME)
            
        end
    end
end