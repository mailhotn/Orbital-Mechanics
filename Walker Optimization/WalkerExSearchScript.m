%% Initialization
clear
PropParams.maxPdop = 1000;
PropParams.timeVec = 0:100:86164;
PropParams.elevMin  = 5;
PropParams.relTol = 1e-6;
PropParams.absTol = 1e-6;
PropParams.datafolder = 'C:\Users\User\Dropbox\Walker Optimization Data';
PropParams.delLat = 5;
primary = earth();

nRepeats = 14;
nDays    = 1;

latList = 30:10:60;

% latEm = 40;
lonEm = 0;
maxSats = 80;
minSats = 40;

dIncList = 10; % Option for optimizing inclination, not fully implemented
% hAList = dIncList;
hAList = 10;
save([PropParams.datafolder '\OptParams.mat']);

%% Massive For Loop

parfor iLat = 1:length(latList)
    latEm = latList(iLat);
    for iDInc = 1:numel(dIncList)
        Orbit = struct();
        [sma, ecc] = CalcRgtSmaApoHeight(latEm+dIncList(iDInc),0, nRepeats, nDays);
        Orbit.sma = sma;
        Orbit.ecc = ecc;
        Orbit.inc = latEm + dIncList(iDInc);
        Orbit.hA = dIncList(iDInc);
        %     Orbit = CalcOptIncRoiA(nRepeats,nDays,latEm,0,5);
        InitCon = InitConElliptical(Orbit.ecc,Orbit.inc,Orbit.sma,latEm,lonEm);
        for nSats = minSats:maxSats
            % Optimize
            try
                % Run Search for T & latGs
                tic
                Arch = struct();
                Arch.nSats = nSats;
                Arch.nRepeats = nRepeats;
                Arch.nDays = nDays;
                ExSol = WalkerExSearch( Arch, Orbit, InitCon, latEm, PropParams);
                optTime = toc;
                % Verbose Output Message
                c = clock;
                disp([newline num2str(c(3)) '/' num2str(c(2)) ' ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6),2)...
                    newline 'Optimization Complete for T = ' num2str(nSats) ...
                    ' and Latitude = ' num2str(latEm)...
                    newline 'Fitness: ' num2str(ExSol.fit)...
                    newline 'Elapsed Time: ' num2str((optTime-mod(optTime,60))/60) ' min '...
                    num2str(mod(optTime,60),3) ' sec'...
                    newline 'Optimal Solution i:T/P/F: '...
                    num2str(ExSol.orbits{ExSol.iOpt}.inc) '°'  ':' num2str(nSats) '/' num2str(ExSol.optNPlanes)...
                    '/' num2str(-ExSol.phaseMat(1,ExSol.iOpt)) ...
                    newline 'Semimajor Axis: ' num2str(ExSol.orbits{ExSol.iOpt}.sma) ' km' ...
                    newline 'Earth Repeats: ' num2str(nRepeats)])
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
end


c = clock;
fileID = fopen(['C:\Users\User\Dropbox\zzzMatlab Notification\OptDone_'...
    num2str(c(3)) '-' num2str(c(2)) '_' num2str(c(4)) ...
    '.' num2str(c(5)) '.txt'],'w'); % create error file for IFTT phone notification
fprintf(fileID,'done');