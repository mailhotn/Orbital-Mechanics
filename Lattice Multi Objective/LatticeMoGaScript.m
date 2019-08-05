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

latList = 60;
maxSats = 100;
minSats = 45;

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
            'PopulationSize',400,...
            'UseParallel', true,...
            'CrossoverFraction',0.4,...
            'ParetoFraction',0.2,...
            'StallGenLimit',200,...
            'TolFun',5e-6,...
            'OutputFcn', @(options,state,flag) ...
            GaMoGenSave(options,state,flag,OptParams,latEm),...
            'PlotFcns',{@gaplotpareto, @gaplotspread,...
            @gaplotscorediversity, @gaplotrankhist, @gaplotstopping});
        
        if exist([OptParams.datafolder '\MoGaTemp_Lat_' num2str(latEm)...
                '.mat'],'file')
            % Load from autosave if there is one
            load([OptParams.datafolder '\MoGaTemp_Lat_' num2str(latEm)...
                '.mat']);
            options = gaoptimset(options,'InitialPopulation',state.Population);
            disp(['Autosave found for Latitude ' num2str(latEm)...
                newline 'Continuing from Generation ' num2str(state.Generation)])
            clear state;            
        end
        
        [x, f, exitflag, output, ~, ~] = ...
            gamultiobj(@(x) LatMoGaFitness(x,latEm,OptParams),...
            7, [], [], [], [], lb, ub, options);
        nPareto = size(x,1);
        x = sortrows(x);
        GaSol = struct();
        for iPar = 1:nPareto
            parSol = LatMoGaFitnessPlus(x(iPar,:),latEm,OptParams);
            GaSol.Cons{iPar} = parSol.Con;
            GaSol.coverage(:,iPar) = parSol.coverage;
            GaSol.maxPdop(:,iPar) = parSol.maxPdop;
            GaSol.intPdop(:,iPar) = parSol.intPdop;
            GaSol.p90(:,iPar) = parSol.p90;
            GaSol.p75(:,iPar) = parSol.p75;
            GaSol.p50(:,iPar) = parSol.p50;
            GaSol.latEm = latEm;
            GaSol.OptParams = OptParams;
        end
        save([OptParams.datafolder '\LatticeGaSol_Lat_' num2str(latEm)...
            '.mat'],'GaSol');
        delete([OptParams.datafolder '\MoGaTemp_Lat_' num2str(latEm)...
                '.mat']);
        
        optTime = toc;
        % Verbose Output Message
        c = clock;
        disp([newline num2str(c(3)) '/' num2str(c(2)) ' ' num2str(c(4)) ':' num2str(c(5)) ':' num2str(c(6),2)...
            newline 'Optimization Complete for Latitude = ' num2str(latEm) ...
            newline 'Elapsed Time: ' num2str((optTime-mod(optTime,60))/60) ' min '...
            num2str(mod(optTime,60),3) ' sec'...
            newline 'Generations: ' num2str(output.generations)...
            newline 'Pareto Individuals: ' num2str(nPareto)])
        
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