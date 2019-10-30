function GaSol = LatticeGa(Arch, latEm, OptParams)
% Run GA for each combination of Ns No in order to find optimal hA, inc,
% phasing
if ~exist([OptParams.datafolder '\LatticeGaSol_Lat_' num2str(latEm)...
        '_nSats_' num2str(Arch.nSats) '.mat'],'file')
    pList = divisors(Arch.nSats);
    
    % initialize Performance Arrays
    maxMat = inf(3,length(pList));
    intMat = inf(3,length(pList));
    p90Mat = inf(3,length(pList));
    p75Mat = inf(3,length(pList));
    p50Mat = inf(3,length(pList));
    covMat = zeros(3,length(pList));
    
    % Initialize Constellation Parameter Arrays
    phaseParams = nan(3,length(pList)); % nC1, nC2, nC3
    archParams  = nan(3,length(pList)); % nPlanes, nAops, nSatsPerAop
    
    for iPlanes = 1:length(pList)
        % Initialize Ga
        Arch.nPlanes = pList(iPlanes);
        Arch.nSatsPerAop = OptParams.nDays; % Constraint!  can be different if nDays > 1
        Arch.nAops = Arch.nSats/Arch.nPlanes/Arch.nSatsPerAop;
        GaOptions = gaoptimset('PopulationSize',50);
        intCon = [1,3,4,5];
        lBounds = [1,... % hA
            5,... % Inc
            1,... % nC1
            1,... % nC2
            1];   % nC3
        uBounds = [numel(OptParams.hAList),... % hA
            30,...                      % Inc
            Arch.nPlanes,...            % nC1
            Arch.nAops,...              % nC2
            Arch.nPlanes];              % nC3
        % Optimize
        Sol = ga(@(x) LatticeGaFitness(x,Arch,latEm,OptParams),5,[],...
            [],[],[],lBounds,uBounds,[],intCon,GaOptions);
        
        % Re-evaluate (Should prob use exit function instead)
        hA = OptParams.hAList(Sol(1));
        inc = latEm + Sol(2);
        Phase.nC1 = Sol(3);
        Phase.nC2 = Sol(4);
        Phase.nC3 = Sol(5);
        % Orbital Parameters
        [sma,ecc] = CalcRgtSmaApoHeight(inc, hA, OptParams.nRepeats, OptParams.nDays);
        Orbit = struct();
        Orbit.sma = sma;
        Orbit.ecc = ecc;
        Orbit.inc = inc;
        Orbit.hA = hA;
        
        % Initial Conditions
        InitCon = InitConElliptical(ecc,inc,sma,latEm,0);
        
        % Create Constellation & Propagate
        LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);
        Prop = Propagator(LC,OptParams.relTol,OptParams.absTol);
        [propTime, propState] = Prop.PropEciJ2(OptParams.timeVec);
        
        % Evaluate PDOP
        [pdop, ~] = TdoaPdopVec(propState,propTime,latEm,0,0,OptParams.elevMin);
        [pdopN, ~] = TdoaPdopVec(propState,propTime,...
            latEm+OptParams.delLat,0,0,OptParams.elevMin);
        [pdopS, ~] = TdoaPdopVec(propState,propTime,...
            latEm-OptParams.delLat,0,0,OptParams.elevMin);
        pdopMat = [pdopN.';pdop.';pdopS.'];
        % Evaluate Performance at all latitudes
        for iLat = 1:3
            coverage = 100 - sum(isnan(pdopMat(iLat,:)))/length(pdopMat(iLat,:))*100;
            if any(~isnan(pdopMat(iLat,:)))
                maxPdop  = max(pdopMat(iLat,~isnan(pdopMat(iLat,:))));
                pdopMat(iLat,pdopMat(iLat,:) > OptParams.maxPdop) = OptParams.maxPdop;
                pdopMat(iLat,isnan(pdopMat(iLat,:))) = OptParams.maxPdop;
                intPdop  = trapz(propTime,pdopMat(iLat,:))/(propTime(end)-propTime(1));
                intMat(iLat,iPlanes) = intPdop;
                maxMat(iLat,iPlanes) = maxPdop;
                covMat(iLat,iPlanes) = coverage;
                p90Mat(iLat,iPlanes) = prctile(pdopMat(iLat,:),90);
                p75Mat(iLat,iPlanes) = prctile(pdopMat(iLat,:),75);
                p50Mat(iLat,iPlanes) = prctile(pdopMat(iLat,:),50);
                phaseParams(:,iPlanes) = [Phase.nC1;Phase.nC2;Phase.nC3];
                archParams(:,iPlanes)  = [Arch.nPlanes;Arch.nAops;Arch.nSatsPerAop];
            end
        end
    end
    [fit,iOpt] = min(intMat(2,:));
    GaSol.nSats = Arch.nSats;
    GaSol.archMat  = archParams;
    GaSol.phaseMat = phaseParams;
    GaSol.Orbit = Orbit;
    GaSol.InitCon = InitCon;
    GaSol.latEm = latEm;
    GaSol.optNPlanes = archParams(1,iOpt);
    GaSol.iOpt = iOpt;
    GaSol.OptParams = OptParams;
    GaSol.coverage = covMat;
    GaSol.maxPdop  = maxMat;
    GaSol.intPdop  = intMat;
    GaSol.p90 = p90Mat;
    GaSol.p75 = p75Mat;
    GaSol.p50 = p50Mat;
    GaSol.fit = fit;
    save([OptParams.datafolder '\LatticeGaSol_Lat_' num2str(latEm)...
        '_nSats_' num2str(Arch.nSats) '.mat'],'GaSol');
else
    load([OptParams.datafolder '\LatticeGaSol_Lat_' num2str(latEm)...
        '_nSats_' num2str(Arch.nSats) '.mat']);
end         