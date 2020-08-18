clear
targetFolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\GA Standard\Previous Runs\Version 6 - definitive';
sourceFolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\GA Standard\Previous Runs\Version 5 - dT 100, Pop 80';
load([sourceFolder '\OptParams.mat'],'OptParams','latList',...
    'minSats','maxSats');
OptParams.minMinDist = 1; % km
save([targetFolder '\OptParams.mat'],'OptParams','latList',...
    'minSats','maxSats');
%%
parfor iLat = 1:length(latList)
    latEm = latList(iLat);
    for nSats = minSats:maxSats
        try
            tic
            [GaSol,nColls] = LatticeGaResColl(latEm,nSats,...
                OptParams,sourceFolder,targetFolder);
            optTime = toc;
            % Verbose Output Message
            c = clock;
            disp([newline num2str(c(3)) '/' num2str(c(2)) ' ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6),2)...
                newline 'Collision Checking Complete for nSats = ' num2str(nSats) ...
                ', Latitude = ' num2str(latEm)...
                newline 'Collisions Resolved: ' num2str(nColls)...
                newline 'Elapsed Time: ' num2str((optTime-mod(optTime,60))/60) ' min '...
                num2str(mod(optTime,60),3) ' sec'])
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
