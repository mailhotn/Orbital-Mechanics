function [GaSol,nColls] = LatticeGaResColl(latEm,nSats,OptParams,sourceFolder,targetFolder)

nColls = 0;
if ~exist([targetFolder '\LatticeGaSol_Lat_' num2str(latEm)...
        '_nSats_' num2str(nSats) '.mat'],'file')
    
    load([sourceFolder '\LatticeGaSol_Lat_' num2str(latEm) '_nSats_' ...
        num2str(nSats) '.mat'],'GaSol');
    
    minDists = zeros(1,size(GaSol.archMat,2));
    
    for iPlanes = 1:size(GaSol.archMat,2)
        Arch.nSats = nSats;
        Arch.nPlanes = GaSol.archMat(1,iPlanes);
        Arch.nAops = GaSol.archMat(2,iPlanes);
        Arch.nSatsPerAop = GaSol.archMat(3,iPlanes);
        Arch.nRepeats = 14;
        Arch.nDays = 1;
        
        Phase.nC1 = GaSol.phaseMat(1,iPlanes);
        Phase.nC2 = GaSol.phaseMat(2,iPlanes);
        Phase.nC3 = GaSol.phaseMat(3,iPlanes);
        
        Orbit = GaSol.orbits{iPlanes};
        InitCon = GaSol.inits{iPlanes};
        
        LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);
        
        minDist = CalcMinDist(LC);
        
        if minDist < OptParams.minMinDist   
            nColls = nColls + 1; % Collision Detected!!!
            GaOptions = gaoptimset('PopulationSize',80);
            intCon = [1,3,4,5];
            lBounds = [1,... % hA
                5,... % Inc
                1,... % nC1
                1,... % nC2
                1];   % nC3
            uBounds = [numel(OptParams.hAList),... % hA
                20,...                      % Inc
                Arch.nPlanes,...            % nC1
                Arch.nAops,...              % nC2
                Arch.nPlanes];              % nC3
            % Optimize
            Sol = ga(@(x) LatticeGaFitness(x,Arch,latEm,OptParams),5,[],...
                [],[],[],lBounds,uBounds,@(x) LatticeGaCollCon(x,Arch,latEm,OptParams)...
                ,intCon,GaOptions);
            
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
            minDist = CalcMinDist(LC);
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
                    GaSol.intPdop(iLat,iPlanes) = intPdop;
                    GaSol.maxPdop(iLat,iPlanes) = maxPdop;
                    GaSol.coverage(iLat,iPlanes) = coverage;
                    GaSol.p90(iLat,iPlanes) = prctile(pdopMat(iLat,:),90);
                    GaSol.p75(iLat,iPlanes) = prctile(pdopMat(iLat,:),75);
                    GaSol.p50(iLat,iPlanes) = prctile(pdopMat(iLat,:),50);
                    GaSol.phaseMat(:,iPlanes) = [Phase.nC1;Phase.nC2;Phase.nC3];
                    GaSol.archMat(:,iPlanes)  = [Arch.nPlanes;Arch.nAops;Arch.nSatsPerAop];
                    GaSol.orbits{iPlanes} = Orbit;
                    GaSol.inits{iPlanes} = InitCon;
                end
            end
        end
        minDists(iPlanes) = minDist;
    end
    [fit,iOpt] = min(GaSol.intPdop(2,:));
    GaSol.iOpt = iOpt;
    GaSol.fit = fit;
else
    load([targetFolder '\LatticeGaSol_Lat_' num2str(latEm)...
        '_nSats_' num2str(nSats) '.mat'],'GaSol');
end
save([targetFolder '\LatticeGaSol_Lat_' num2str(latEm)...
    '_nSats_' num2str(nSats) '.mat'],'GaSol');