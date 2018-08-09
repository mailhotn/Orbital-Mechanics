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
primary = earth();
nOrbits = 6;

maxSats = 80;
minSats = 20;

delInc = 25;

incList = [5,10,20,25];

% save([PropParams.datafolder '\OptParams.mat']);
%% Massive For Loop

parfor iLat = 1:length(latList)
    GroundS = struct();
    GroundS.lat = latList(iLat);
    GroundS.lon = 0;
    ecc = 0;
%     [inc, sma, raan0] = VirtualSatOe(ecc,Rgt,GroundS,nOrbits);
    inc = min([90,latList(iLat)+delInc]);
    sma = CalcRgtSma(ecc,inc,Rgt.jRepeats,Rgt.kDays);
    
    n = sqrt(primary.mu/sma^3);
    raanDot = -1.5*(primary.J2*sqrt(primary.mu)*primary.Re^2)/(sma^3.5*(1-ecc^2)^2)*cosd(inc);
    aol0 = 0;
    t0 = 0;
    aol1 = asin(sind(GroundS.lat)/sind(inc));
    t1 = (aol1-aol0)/n + t0;
    
    raan0 = GroundS.lon + (primary.we - raanDot*180/pi)*(t1-t0) - ...
        atan2d(cosd(inc),(tan(aol1 - aol0))^-1);
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
            fileID = fopen(['C:\Users\User\Dropbox\zzzMatlab Errors\error_'...
                num2str(c(3)) '-' num2str(c(2)) '_' num2str(c(4)) ...
                '.' num2str(c(5)) '.txt'],'w'); % create error file for IFTT phone notification
            fprintf(fileID,ME.message);
            rethrow(ME)
            
        end
    end
end

c = clock;
fileID = fopen(['C:\Users\User\Dropbox\zzzMatlab Notification\OptDone_'...
    num2str(c(3)) '-' num2str(c(2)) '_' num2str(c(4)) ...
    '.' num2str(c(5)) '.txt'],'w'); % create error file for IFTT phone notification
fprintf(fileID,'done');