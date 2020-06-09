clear
datafolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\GA Standard\Previous Runs\Version 2 Copy for Inc restoration';
load([datafolder '\OptParams.mat']);
%% Ga Test
load([datafolder '\LatticeGaSol_Lat_30_nSats_60.mat'])
nCons = size(GaSol.archMat,2);

pList = divisors(Arch.nSats);

% initialize Performance Arrays
maxMat = inf(3,length(pList));
intMat = inf(3,length(pList));
p90Mat = inf(3,length(pList));
p75Mat = inf(3,length(pList));
p50Mat = inf(3,length(pList));
covMat = zeros(3,length(pList));

error = inf(1,length(pList));


% Initialize Constellation Parameter Arrays
phaseParams = nan(3,length(pList)); % nC1, nC2, nC3
archParams  = nan(3,length(pList)); % nPlanes, nAops, nSatsPerAop
orbits = cell(1,length(pList));
inits  = cell(1,length(pList));

tic
for iPlanes = 1:nCons
    Arch = struct();
    Arch.nSats = GaSol.nSats;
    Arch.nPlanes = GaSol.archMat(1,iPlanes);
    Arch.nAops = GaSol.archMat(2,iPlanes);
    Arch.nSatsPerAop = GaSol.archMat(3,iPlanes);
    Phase = struct();
    Phase.nC1 = GaSol.phaseMat(1,iPlanes);
    Phase.nC2 = GaSol.phaseMat(2,iPlanes);
    Phase.nC3 = GaSol.phaseMat(3,iPlanes);

    % Initialize Ga
    GaOptions = gaoptimset('PopulationSize',10,'UseParallel',true); %,...
%         'PlotFcn',{'gaplotdistance','gaplotrange'});
    intCon = 1;
    lBounds = [1,... % hA
        5];   % Inc
    uBounds = [numel(OptParams.hAList),... % hA
        30];                        % Inc
    % Optimize
    Sol = ga(@(x) LatticeIncRestGaFitness(x,Arch,Phase,GaSol.latEm,OptParams),2,[],...
        [],[],[],lBounds,uBounds,[],intCon,GaOptions);
    
    hA = OptParams.hAList(Sol(1));
    inc = GaSol.latEm + Sol(2);
    % Orbital Parameters
    [sma,ecc] = CalcRgtSmaApoHeight(inc, hA, OptParams.nRepeats, OptParams.nDays);
    Orbit = struct();
    Orbit.sma = sma;
    Orbit.ecc = ecc;
    Orbit.inc = inc;
    Orbit.hA = hA;
    
    % Initial Conditions
    InitCon = InitConElliptical(ecc,inc,sma,GaSol.latEm,0);
    
    % Create Constellation & Propagate
    LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);
    Prop = Propagator(LC,OptParams.relTol,OptParams.absTol);
    [propTime, propState] = Prop.PropEciJ2(OptParams.timeVec);
    
    % Evaluate PDOP
    [pdop, ~] = TdoaPdopVec(propState,propTime,GaSol.latEm,0,0,OptParams.elevMin);
    [pdopN, ~] = TdoaPdopVec(propState,propTime,...
        GaSol.latEm+OptParams.delLat,0,0,OptParams.elevMin);
    [pdopS, ~] = TdoaPdopVec(propState,propTime,...
        GaSol.latEm-OptParams.delLat,0,0,OptParams.elevMin);
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
            orbits{iPlanes} = Orbit;
            inits{iPlanes} = InitCon;
        end
    end
    error(iPlanes) = norm(intMat(:,iPlanes)-GaSol.intPdop(:,iPlanes));
end
toc
error


%% fminbndtest

load([datafolder '\LatticeGaSol_Lat_30_nSats_70.mat'])
nCons = size(GaSol.archMat,2);
tic
pList = divisors(GaSol.nSats);

% initialize Performance Arrays
maxMat = inf(3,length(pList));
intMat = inf(3,length(pList));
p90Mat = inf(3,length(pList));
p75Mat = inf(3,length(pList));
p50Mat = inf(3,length(pList));
covMat = zeros(3,length(pList));

error = inf(1,length(pList));

% Initialize Constellation Parameter Arrays
phaseParams = nan(3,length(pList)); % nC1, nC2, nC3
archParams  = nan(3,length(pList)); % nPlanes, nAops, nSatsPerAop
orbits = cell(1,length(pList));
inits  = cell(1,length(pList));
for iPlanes = 1:nCons
    Arch = struct();
    Arch.nSats = GaSol.nSats;
    Arch.nPlanes = GaSol.archMat(1,iPlanes);
    Arch.nAops = GaSol.archMat(2,iPlanes);
    Arch.nSatsPerAop = GaSol.archMat(3,iPlanes);
    Phase = struct();
    Phase.nC1 = GaSol.phaseMat(1,iPlanes);
    Phase.nC2 = GaSol.phaseMat(2,iPlanes);
    Phase.nC3 = GaSol.phaseMat(3,iPlanes);
    

    fit = inf;
    for iHA = 1:numel(OptParams.hAList)
        hA = OptParams.hAList(iHA);
        [inc,fval] = fminbnd(@(x) LatticeIncRestFminFitness(x,Arch,Phase,hA,...
            GaSol.latEm, GaSol.intPdop(:,iPlanes), OptParams),8,20);
        inc = inc + GaSol.latEm;
        if fval<fit
            fit = fval;
            % Orbital Parameters
            [sma,ecc] = CalcRgtSmaApoHeight(inc, hA, OptParams.nRepeats, OptParams.nDays);
            Orbit = struct();
            Orbit.sma = sma;
            Orbit.ecc = ecc;
            Orbit.inc = inc;
            Orbit.hA = hA;
            
            % Initial Conditions
            InitCon = InitConElliptical(ecc,inc,sma,GaSol.latEm,0);
            
            % Create Constellation & Propagate
            LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);
            Prop = Propagator(LC,OptParams.relTol,OptParams.absTol);
            [propTime, propState] = Prop.PropEciJ2(OptParams.timeVec);
            
            % Evaluate PDOP
            [pdop, ~] = TdoaPdopVec(propState,propTime,GaSol.latEm,0,0,OptParams.elevMin);
            [pdopN, ~] = TdoaPdopVec(propState,propTime,...
                GaSol.latEm+OptParams.delLat,0,0,OptParams.elevMin);
            [pdopS, ~] = TdoaPdopVec(propState,propTime,...
                GaSol.latEm-OptParams.delLat,0,0,OptParams.elevMin);
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
                    orbits{iPlanes} = Orbit;
                    inits{iPlanes} = InitCon;
                end
            end
            error(iPlanes) = norm(intMat(:,iPlanes)-GaSol.intPdop(:,iPlanes));
        end
    end
end
toc
error

%% pattern test

load([datafolder '\LatticeGaSol_Lat_30_nSats_60.mat'])
nCons = size(GaSol.archMat,2);
tic
pList = divisors(Arch.nSats);

% initialize Performance Arrays
maxMat = inf(3,length(pList));
intMat = inf(3,length(pList));
p90Mat = inf(3,length(pList));
p75Mat = inf(3,length(pList));
p50Mat = inf(3,length(pList));
covMat = zeros(3,length(pList));

error = inf(1,length(pList));

% Initialize Constellation Parameter Arrays
phaseParams = nan(3,length(pList)); % nC1, nC2, nC3
archParams  = nan(3,length(pList)); % nPlanes, nAops, nSatsPerAop
orbits = cell(1,length(pList));
inits  = cell(1,length(pList));
for iPlanes = 1:nCons
    Arch = struct();
    Arch.nSats = GaSol.nSats;
    Arch.nPlanes = GaSol.archMat(1,iPlanes);
    Arch.nAops = GaSol.archMat(2,iPlanes);
    Arch.nSatsPerAop = GaSol.archMat(3,iPlanes);
    Phase = struct();
    Phase.nC1 = GaSol.phaseMat(1,iPlanes);
    Phase.nC2 = GaSol.phaseMat(2,iPlanes);
    Phase.nC3 = GaSol.phaseMat(3,iPlanes);
    

    fit = inf;
    for iHA = 1:numel(OptParams.hAList)
        hA = OptParams.hAList(iHA);
        [inc,fval] = patternsearch(@(x) LatticeIncRestFminFitness(x,Arch,Phase,hA,...
            GaSol.latEm,OptParams),10,[],[],[],[],5,30);
        inc = inc + GaSol.latEm;
        if fval<fit
            fit = fval;
            % Orbital Parameters
            [sma,ecc] = CalcRgtSmaApoHeight(inc, hA, OptParams.nRepeats, OptParams.nDays);
            Orbit = struct();
            Orbit.sma = sma;
            Orbit.ecc = ecc;
            Orbit.inc = inc;
            Orbit.hA = hA;
            
            % Initial Conditions
            InitCon = InitConElliptical(ecc,inc,sma,GaSol.latEm,0);
            
            % Create Constellation & Propagate
            LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);
            Prop = Propagator(LC,OptParams.relTol,OptParams.absTol);
            [propTime, propState] = Prop.PropEciJ2(OptParams.timeVec);
            
            % Evaluate PDOP
            [pdop, ~] = TdoaPdopVec(propState,propTime,GaSol.latEm,0,0,OptParams.elevMin);
            [pdopN, ~] = TdoaPdopVec(propState,propTime,...
                GaSol.latEm+OptParams.delLat,0,0,OptParams.elevMin);
            [pdopS, ~] = TdoaPdopVec(propState,propTime,...
                GaSol.latEm-OptParams.delLat,0,0,OptParams.elevMin);
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
                    orbits{iPlanes} = Orbit;
                    inits{iPlanes} = InitCon;
                end
            end
            error(iPlanes) = norm(intMat(:,iPlanes)-GaSol.intPdop(:,iPlanes));
        end
    end
end
toc
error