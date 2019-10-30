%% Initialize Optimization Parameters
clear
OptParams.maxPdop = 1000;
OptParams.timeVec = 0:10:86164;
OptParams.elevMin = 5;
OptParams.relTol  = 1e-6;
OptParams.absTol  = 1e-6;
OptParams.datafolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\GA Standard';
OptParams.nRepeats = 14;
OptParams.nDays    = 1;
OptParams.hAList = [0, 900, 1000];
OptParams.delLat = 5;

latList = 60;
maxSats = 80;
minSats = 40;

save([OptParams.datafolder '\OptParams.mat']);
%% Run Genetic Algorithm
for iLat = 1:length(latList)
    latEm = latList(iLat);
    parfor nSats = minSats:maxSats
        % Optimize
        try
            tic
            Arch = struct();
            Arch.nSats = nSats;
            Arch.nRepeats = OptParams.nRepeats;
            Arch.nDays = OptParams.nDays;
            
            GaSol = LatticeGa(Arch,latEm,OptParams);
            
            optTime = toc;
            % Verbose Output Message
            c = clock;
            disp([newline num2str(c(3)) '/' num2str(c(2)) ' ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6),2)...
                newline 'Optimization Complete for nSats = ' num2str(nSats) ...
                ' and Latitude = ' num2str(latEm)...
                newline 'Fitness: ' num2str(GaSol.fit)...
                newline 'Elapsed Time: ' num2str((optTime-mod(optTime,60))/60) ' min '...
                num2str(mod(optTime,60),3) ' sec'...
                newline 'Optimal Solution nSats/nPlanes/nC1/nC2/nC3: '...
                num2str(nSats) '/' num2str(GaSol.archMat(1,GaSol.iOpt))...
                '/' num2str(GaSol.phaseMat(1,GaSol.iOpt)) ...
                '/' num2str(GaSol.phaseMat(2,GaSol.iOpt)) ...
                '/' num2str(GaSol.phaseMat(3,GaSol.iOpt)) ...
                newline 'Inclination: ' num2str(GaSol.Orbit.inc) '°'...
                newline 'Eccentricity: ' num2str(GaSol.Orbit.ecc)...
                newline 'Semimajor Axis: ' num2str(GaSol.Orbit.sma) ' km' ...
                newline 'Earth Repeats: ' num2str(OptParams.nRepeats)])
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