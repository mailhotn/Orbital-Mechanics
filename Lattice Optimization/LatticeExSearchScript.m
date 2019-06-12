%% Initialize Optimization Parameters
clear
PropParams.maxPdop = 1000;
PropParams.timeVec = 0:10:86164;
PropParams.elevMin = 5;
PropParams.relTol  = 1e-6;
PropParams.absTol  = 1e-6;
PropParams.datafolder = 'C:\Users\User\Dropbox\Lattice Optimization Data';
primary = earth();

nRepeats = 14;
nDays    = 1;

latList = 10:10:80;

% latEm = 40;
lonEm = 0;
maxSats = 120;
minSats = 15;

hAList = [0, 900, 1000];

save([PropParams.datafolder '\OptParams.mat']);

%% Perform Search
parfor iLat = 1:length(latList)

    latEm = latList(iLat)
    for iHA = 1:length(hAList)
        Orbit = struct();
        [sma, ecc, inc] = RgtSunSynElements(hAList(iHA), nRepeats, nDays);
        Orbit.sma = sma;
        Orbit.ecc = ecc;
        Orbit.inc = inc;
        Orbit.hA = hAList(iHA);
        
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

                
                ExSol = LatticeExSearch(Arch,Orbit,InitCon,latEm,PropParams);

                optTime = toc;
                % Verbose Output Message
                c = clock;
                disp([newline num2str(c(3)) '/' num2str(c(2)) ' ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6),2)...
                    newline 'Optimization Complete for nSats = ' num2str(nSats) ...
                    ' and Latitude = ' num2str(latEm)...
                    newline 'Fitness: ' num2str(ExSol.fit)...
                    newline 'Elapsed Time: ' num2str((optTime-mod(optTime,60))/60) ' min '...
                    num2str(mod(optTime,60),3) ' sec'...
                    newline 'Optimal Solution nSats/nPlanes/nC1/nC2/nC3: '...
                    num2str(nSats) '/' num2str(ExSol.archMat(1,ExSol.iOpt))...
                    '/' num2str(ExSol.phaseMat(1,ExSol.iOpt)) ...
                    '/' num2str(ExSol.phaseMat(2,ExSol.iOpt)) ...
                    '/' num2str(ExSol.phaseMat(3,ExSol.iOpt)) ...
                    newline 'Inclination: ' num2str(Orbit.inc) '°'...
                    newline 'Eccentricity: ' num2str(Orbit.ecc)...
                    newline 'Semimajor Axis: ' num2str(Orbit.sma) ' km' ...
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