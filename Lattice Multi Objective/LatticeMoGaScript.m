%% Initialize Optimization Parameters
clear
OptParams.maxPdop = 1000;
OptParams.timeVec = 0:10:86164;
OptParams.elevMin = 5;
OptParams.relTol  = 1e-6;
OptParams.absTol  = 1e-6;
OptParams.datafolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Multi Objective';
OptParams.nRepeats = 14;
OptParams.nDays    = 1;
OptParams.delLat = 5;

latList = 40;
maxSats = 80;
minSats = 40;

save([OptParams.datafolder '\OptParams.mat']);
%% Run Genetic Algorithm
for iLat = 1:length(latList)
    latEm = latList(iLat);
    % Optimize
    try
        tic
        lb = [minSats, -1e6, -1e6, -1e6, -1e6, latEm, -1e6];
        ub = [maxSats, 1e6, 1e6, 1e6, 1e6, min([latEm+30, 180-latEm]), 1e6];
        bounds = [lb; ub];
        options = gaoptimset('CreationFcn', @LatIntCreation,...
            'MutationFcn', @LatIntMutation,...
            'CrossoverFcn', @LatIntCrossover,...
            'PopInitRange', bounds,...
            'PopulationSize',100,...
            'UseParallel', false,...
            'PlotFcns',{@gaplotpareto, @gaplotparetodistance,...
            @gaplotscorediversity, @gaplotgenealogy, @gaplotrankhist});
        
        [x, f, exitflag, output, population, score] = ...
            gamultiobj(@(x) LatMoGaFitness(x,latEm,OptParams),...
            7, [], [], [], [], lb, ub, options);
        
        optTime = toc;
        
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


c = clock;
fileID = fopen(['C:\Users\User\Dropbox\zzzMatlab Notification\OptDone_'...
    num2str(c(3)) '-' num2str(c(2)) '_' num2str(c(4)) ...
    '.' num2str(c(5)) '.txt'],'w'); % create error file for IFTT phone notification
fprintf(fileID,'done');