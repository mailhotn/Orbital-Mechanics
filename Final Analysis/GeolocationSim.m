%% Define Constellation & RoI Position
clear
latRoi = 30;
lonRoi = 35;
delLatRoi = 10;
delLonRoi = 10;
elevMin = 5; % minimum elevation for visibility [deg]
load('C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\LatticeDef v1\LatticeExSol_Lat_30_nSats_70_hA_900.mat');

Arch.nSats = 70;
Arch.nPlanes = 7;
iCon = find(ExSol.archMat(1,:)==Arch.nPlanes);
Arch.nSatsPerAop = ExSol.archMat(3,iCon);
Arch.nAops = ExSol.archMat(2,iCon);

Phase.nC1 = ExSol.phaseMat(1,iCon);
Phase.nC2 = ExSol.phaseMat(2,iCon);
Phase.nC3 = ExSol.phaseMat(3,iCon);

Orbit = ExSol.orbits{iCon};

InitCon.raan1 = lonRoi-atand(cosd(Orbit.inc)*tand(latRoi));
InitCon.aop1 = asind(sind(latRoi)/sind(Orbit.inc)) - 180;
InitCon.M1 = 180;
% InitCon = ExSol.inits{iCon};


LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);

% figure(1)
% PlotGroundTrack(LC);
figure(2)
PlotPdopMap(LC,latRoi,lonRoi,5,true,true,95);
hold on
plot([lonRoi - 40,lonRoi + 40],[latRoi + delLatRoi, latRoi + delLatRoi],'-w',...
    [lonRoi - 40,lonRoi + 40],[latRoi - delLatRoi, latRoi - delLatRoi],'-w','linewidth',2);
hold off
%% Orbital Simulation
time = 0:600:86400; % sec
Prop = Propagator(LC,1e-8,1e-8);
[~, state] = Prop.PropEciJ2(time);

%% Estimation
nPoints = 1000;
nExps = 100;

latList = nan(nPoints,1);
lonList = nan(nPoints,1);

rUnConAbsErr = nan(nPoints,length(time));
rLinConAbsErr = nan(nPoints,length(time));
r2ndConAbsErr = nan(nPoints,length(time));

rUnConRmsErr = nan(nPoints,length(time));
rLinConRmsErr = nan(nPoints,length(time));
r2ndConRmsErr = nan(nPoints,length(time));

pdopMat = nan(nPoints,length(time));
unCpdopMat = nan(nPoints,length(time));
varTime = (10/3e5)^2; % sec

for iSim = 1:nPoints
    latEm = latRoi + (rand-0.5)*2*delLatRoi;
    lonEm = lonRoi + (rand-0.5)*2*delLonRoi;
    latList(iSim) = latEm;
    lonList(iSim) = lonEm;
    
    [pdopMat(iSim,:), ~, unCpdopMat(iSim,:)] = CumTdoaPdopVec(state,time,latEm,lonEm,0,5);
    rEm = lla2ecef([latEm,lonEm,0]).'/1000;
    % [pdop,nSis] = TdoaPdopVec(state,time,latEm,lonEm,0,elevMin);
    rUnCon = cell(nExps,1);
    
    dirInit = 2*pi*rand;
    latInit = latEm + sin(dirInit)*5;
    lonInit = lonEm + cos(dirInit)*5;
    rInit = lla2ecef([latInit,lonInit,0]).'/1000;
        
    for iExp = 1:nExps
        
        rUnCon{iExp} = nan(3,length(time));
        rUnCon{iExp}(:,1) = rInit; % Initial guess 30°N, 30°E
        rLinCon = rUnCon;
        rNLinCon = rUnCon;
        % Initial Abs Error
        rUnConAbsErr(iSim,1) = norm(rEm - rInit);
        rLinConAbsErr(iSim,1) = norm(rEm - rInit);
        r2ndConAbsErr(iSim,1) = norm(rEm - rInit);
    end
    zMeas = cell(nExps,length(time));
    % Create Measurements
    for iTime = 1:length(time)
        % Derive measurement matrix from satellite positions
        gmst = time(iTime)*LC.primary.we;
        stateEcef = eci2ecef(reshape(state(iTime,:),6,LC.nSats),gmst);
        stateInSight = SatsInSight(stateEcef,rEm,elevMin,norm(rEm));
        rSats = stateInSight(1:3,:);
        nMeas = size(rSats,2)-1;
        % Create Measurement & Linearized matrix
        z = tdoaMeas(rSats,rEm);
        noiseCov = varTime*(eye(nMeas) + ones(nMeas));
        for iExp = 1:nExps
            zMeas{iExp,iTime} = z + mvnrnd(zeros(nMeas,1),noiseCov).'; % add noise
        end
    end
    
    %% Unconstrained
    pEst = cell(nExps,1);
    for iExp = 1:nExps
        pEst{iExp} = 1000*eye(3);
    end
    stdVec = nan(1,length(iTime));
    for iTime = 2:length(time)
        % Derive measurement matrix from satellite positions
        gmst = time(iTime)*LC.primary.we;
        stateEcef = eci2ecef(reshape(state(iTime,:),6,LC.nSats),gmst);
        stateInSight = SatsInSight(stateEcef,rEm,elevMin,norm(rEm));
        rSats = stateInSight(1:3,:);
        nMeas = size(rSats,2)-1;
        noiseCov = varTime*(eye(nMeas) + ones(nMeas));
        
        rSamples = nan(3,nExps);
        for iExp = 1:nExps
            % Create Linearized measurement matrix
            [zEst,hMeas] = tdoaMeasMat(rSats,rUnCon{iExp}(:,iTime-1));
            % EKF
            kGain = (pEst{iExp}*hMeas.')/(hMeas*pEst{iExp}*hMeas.' + noiseCov);
            rUnCon{iExp}(:,iTime) = rUnCon{iExp}(:,iTime-1) + ...
                kGain*(zMeas{iExp,iTime} - zEst);
            pEst{iExp} = (eye(3) - kGain*hMeas)*pEst{iExp}*(eye(3) - kGain*hMeas).' + ...
                kGain*noiseCov*kGain.';
            
            % collect error samples
            rSamples(:,iExp) = rEm - rUnCon{iExp}(:,iTime);
            
        end
        % calculate estimate of estimation error covariance matrix
        rMeanExp = mean(rSamples,2);
        samplePEst = zeros(3);
        for iExp = 1:nExps
%             samplePEst = samplePEst + ...
%                 (rSamples(:,iExp) - rMeanExp)*(rSamples(:,iExp) - rMeanExp).';
            samplePEst = samplePEst + pEst{iExp};
        end
        samplePEst = samplePEst/(nExps-1);
        
        % average errors across experiments
        rUnConAbsErr(iSim,iTime) = mean(vecnorm(rSamples));
        rUnConRmsErr(iSim,iTime) = sqrt(trace(samplePEst));
    end
    
    %% Linearized Constraint
    pEst = cell(nExps,1);
    for iExp = 1:nExps
        pEst{iExp} = 1000*eye(3);
    end
    for iTime = 2:length(time)
        % Derive measurement matrix from satellite positions
        gmst = time(iTime)*LC.primary.we;
        stateEcef = eci2ecef(reshape(state(iTime,:),6,LC.nSats),gmst);
        stateInSight = SatsInSight(stateEcef,rEm,elevMin,norm(rEm));
        rSats = stateInSight(1:3,:);
        nMeas = size(rSats,2)-1;
        noiseCov = varTime*(eye(nMeas) + ones(nMeas));
        for iExp = 1:nExps
            % Create Measurement & Linearized matrix
            [zEst,hMeas] = tdoaMeasMat(rSats,rLinCon{iExp}(:,iTime-1));
            % EKF
            kGain = (pEst{iExp}*hMeas.')/(hMeas*pEst{iExp}*hMeas.' + noiseCov);
            rLinCon{iExp}(:,iTime) = rLinCon{iExp}(:,iTime-1) + ...
                kGain*(zMeas{iExp,iTime} - zEst);
            pEst{iExp} = (eye(3) - kGain*hMeas)*pEst{iExp}*(eye(3) - kGain*hMeas).' + ...
                kGain*noiseCov*kGain.';
            % Apply Constraint
            [rLinCon{iExp}(:,iTime), pEst2] = LinConKF(rLinCon{iExp}(:,iTime)...
                ,rLinCon{iExp}(:,iTime-1),pEst{iExp},pEst{iExp});
            
            
            rSamples(:,iExp) = rEm - rLinCon{iExp}(:,iTime);
            
        end
        rMeanExp = mean(rSamples,2);
        samplePEst = zeros(3);
        for iExp = 1:nExps
            samplePEst = samplePEst + ...
                (rSamples(:,iExp) - rMeanExp)*(rSamples(:,iExp) - rMeanExp).';
        end
        samplePEst = samplePEst/(nExps-1);
        % average errors across experiments
        rLinConAbsErr(iSim,iTime) = mean(vecnorm(rSamples));
        rLinConRmsErr(iSim,iTime) = sqrt(trace(samplePEst));
    end

    
    %% 2nd Order Non-linear Constraint
    pEst = cell(nExps,1);
    for iExp = 1:nExps
        pEst{iExp} = 1000*eye(3);
    end
    for iTime = 2:length(time)
        % Derive measurement matrix from satellite positions
        gmst = time(iTime)*LC.primary.we;
        stateEcef = eci2ecef(reshape(state(iTime,:),6,LC.nSats),gmst);
        stateInSight = SatsInSight(stateEcef,rEm,elevMin,norm(rEm));
        rSats = stateInSight(1:3,:);
        nMeas = size(rSats,2)-1;
        noiseCov = varTime*(eye(nMeas) + ones(nMeas));
        for iExp = 1:nExps
            % Create Measurement & Linearized matrix
            [zEst,hMeas] = tdoaMeasMat(rSats,rNLinCon{iExp}(:,iTime-1));
            % EKF
            kGain = (pEst{iExp}*hMeas.')/(hMeas*pEst{iExp}*hMeas.' + noiseCov);
            rNLinCon{iExp}(:,iTime) = rNLinCon{iExp}(:,iTime-1) + ...
                kGain*(zMeas{iExp,iTime} - zEst);
            pEst{iExp} = (eye(3) - kGain*hMeas)*pEst{iExp}*(eye(3) - kGain*hMeas).' + ...
                kGain*noiseCov*kGain.';
            % Apply constraint
            rNLinCon{iExp}(:,iTime) = NonLinConKF(rNLinCon{iExp}(:,iTime),inv(pEst{iExp}));
           
            
            rSamples(:,iExp) = rEm - rNLinCon{iExp}(:,iTime);
            
        end
        rMeanExp = mean(rSamples,2);
        samplePEst = zeros(3);
        for iExp = 1:nExps
            samplePEst = samplePEst + ...
                (rSamples(:,iExp) - rMeanExp)*(rSamples(:,iExp) - rMeanExp).';
        end
        samplePEst = samplePEst/(nExps-1);
        % average errors across experiments
        r2ndConAbsErr(iSim,iTime) = mean(vecnorm(rSamples));
        r2ndConRmsErr(iSim,iTime) = sqrt(trace(samplePEst));
    end
end
%% Plot Stuff
% AbsError Plot
figure(1)
semilogy(time/3600,mean(rUnConAbsErr,1),'-',...
    time/3600,mean(r2ndConAbsErr,1),'-.',...
    'linewidth',1.5)
xlim([0,24])
xticks(0:2:24)
xlabel('$\rm{Time \left[hr\right]}$','interpreter','latex','fontsize',14)
ylabel('$\| r-\hat{r}\| \left[\rm{km}\right]$','interpreter','latex','fontsize',14)
legend('EKF','CEKF')
grid on

% RMS Error Plot
figure(2)
semilogy(time/3600,mean(rUnConRmsErr,1),'-',...
    time/3600,mean(r2ndConRmsErr,1),'-.',...
    time/3600,mean(unCpdopMat,1)*sqrt(varTime)*3e5,':',...
    time/3600,mean(pdopMat,1)*sqrt(varTime)*3e5,':',...
    'linewidth',2)
xlim([0,24])
xticks(0:2:24)
xlabel('$\rm{Time \left[hr\right]}$','interpreter','latex','fontsize',14)
ylabel('$\textrm{RMSE}\left(\hat{r}\right) \left[\rm{km}\right]$','interpreter','latex','fontsize',14)
legend('EKF','CEKF',...
    'Efficient RMSE - unconstrained','Efficient RMSE - constrained')
grid on

% Final Error Plot
figure(3)
semilogy(1:min([nPoints,100]),rUnConAbsErr(1:min([nPoints,100]),end),'o',...
    1:min([nPoints,100]),rLinConAbsErr(1:min([nPoints,100]),end),'s',...
    1:min([nPoints,100]),r2ndConAbsErr(1:min([nPoints,100]),end),'^')
xlabel('$\textrm{Point \#}$','interpreter','latex','fontsize',14)
ylabel('$\textrm{Final } \| r-\hat{r}\| \left[\rm{km}\right]$','interpreter','latex','fontsize',14)
legend('EKF','CEKF with linearized constraint','CEKF with nonlinear constraint')
ylim([0.7,100])
xlim([1,100])
xticks([1,10:10:100])
grid on

% Final Error Latitude Plot
figure(4)
semilogy(latList,rUnConAbsErr(:,end),'o',...
    latList,rLinConAbsErr(:,end),'s',...
    latList,r2ndConAbsErr(:,end),'^')
xlabel('$\textrm{Latitude}$','interpreter','latex','fontsize',14)
ylabel('$\textrm{Final } \| r-\breve{r}\| \left[\rm{km}\right]$','interpreter','latex','fontsize',14)


% Final Error Longitude Plot
figure(5)
semilogy(lonList,rUnConAbsErr(:,end),'o',...
    lonList,rLinConAbsErr(:,end),'s',...
    lonList,r2ndConAbsErr(:,end),'^')
xlabel('$\textrm{Longitude}$','interpreter','latex','fontsize',14)
ylabel('$\textrm{Final } \| r-\breve{r}\| \left[\rm{km}\right]$','interpreter','latex','fontsize',14)


% % Estimate Map
% unConCo = ecef2lla(rUnCon.'*1000);
% linConCo = ecef2lla(rLinCon.'*1000);
% nLinConCo = ecef2lla(rNLinCon.'*1000);
% figure(2)
% clf
% geoshow('landareas.shp')
% hold on
% plot(unConCo(:,2),unConCo(:,1),...
%     linConCo(:,2),linConCo(:,1),...
%     nLinConCo(:,2),nLinConCo(:,1),'--',...
%     lonEm,latEm,'r*',...
%     'linewidth',1.5)
% axis equal
% xlim([17,37])
% ylim([27,37])
% yticks(25:5:45)
% yticklabels({'25°N','30°N','35°N','40°N','45°N'})
% xticks(15:5:35)
% xticklabels({'15°E','20°E','25°E','30°E','35°E'})
% legend('EKF','Linearized Constraint','Nonlinear Constraint')
% grid on
% hold off

% % Constellation Plot
% figure(3)
% PlotGroundTrack(state,time(1))
%
% % Lat Lon Error
% figure(4)
% semilogy(time/3600,abs(latEm - unConCo(:,1)),...
%     time/3600,abs(latEm - linConCo(:,1)),...
%     time/3600,abs(latEm - nLinConCo(:,1)),...
%     'linewidth',1.5)
% xlim([0,24])
% xticks(0:2:24)
% xlabel('Time [hr]')
% ylabel('$\left|\phi_{em}-\hat{\phi}_{em}\right|$','interpreter','latex')
% legend('EKF','Linearized Constraint','Nonlinear Constraint')
% grid on
%
% % Lat Lon Error
% figure(5)
% semilogy(time/3600,abs(lonEm - unConCo(:,2)),...
%     time/3600,abs(lonEm - linConCo(:,2)),...
%     time/3600,abs(lonEm - nLinConCo(:,2)),...
%     'linewidth',1.5)
% xlim([0,24])
% xticks(0:2:24)
% xlabel('Time [hr]')
% ylabel('$\left|\lambda_{em}-\hat{\lambda}_{em}\right|$','interpreter','latex')
% legend('EKF','Linearized Constraint','Nonlinear Constraint')
% grid on
%
% % Lat Lon Error
% figure(6)
% plot(time/3600,(unConCo(:,3)),...
%     time/3600,(linConCo(:,3)),...
%     time/3600,(nLinConCo(:,3)),...
%     'linewidth',1.5)
% xlim([0,24])
% xticks(0:2:24)
% xlabel('Time [hr]')
% ylabel('$\hat{h}_{em} \left[km\right]$','interpreter','latex')
% legend('EKF','Linearized Constraint','Nonlinear Constraint')
% grid on
%% Functions
function [linConEst, linConP] = LinConKF(unConEst, prevConEst, invW, unConP)
e2 = 0.0818191908426215^2; % eccentricity squared (WGS84)
rE = 6378.137; % Earth Equatorial Radius
conMat = diag([1,1,1/(1-e2)]);
D = prevConEst.'*conMat;
d = rE^2;
invDWD = (D*invW*D.')^-1;
linConEst = unConEst - invW*D.'*invDWD*(D*unConEst - d);
linConP = unConP - invW*D.'*invDWD*D*unConP - unConP*D.'*invDWD*D*invW + ...
    invW*D.'*invDWD*D*unConP*D.'*invDWD*D*invW;
end

function nLinConEst = NonLinConKF(unConEst, W)
e2 = 0.0818191908426215^2; % eccentricity squared (WGS84)
rE = 6378.137; % Earth Equatorial Radius
M = diag([1,1,1/(1-e2)]);
m0 = -rE^2;
C = chol(W);
L = diag([1,1,1/sqrt(1-e2)]);
[~,S,V] = svd(L*inv(C));
eVec = V.'*C*unConEst;
sVec = diag(S);

% Newton-Raphson
lam = 0;
dLam = inf;
lamTol = 1e-10;
iter = 0;
maxIter = 10;
while norm(dLam) > lamTol && iter<maxIter
    q = sum(eVec.^2.*sVec.^2./(1+lam*sVec.^2).^2) + m0;
    dq = -2*sum(eVec.^2.*sVec.^4./(1+lam*sVec.^2).^3);
    dLam = -q/dq;
    lam = lam + dLam;
    iter = iter + 1;
end
nLinConEst = (W + lam*M)\(W*unConEst);
end

function [zMeas,hMat] = tdoaMeas(rSats,rEm)
c = 3e5; % speed of light km/s
dMat = rEm - rSats;
dMag = vecnorm(dMat);
zMeas = 1/c*(dMag(2:end) - dMag(1)).'; % True measurement
dMatNorm = dMat./dMag;
hMat = 1/c*(dMatNorm(:,2:end) - dMatNorm(:,1)).';
end

function [zEst,hMatEst] = tdoaMeasMat(rSats,rEst)
c = 3e5; % speed of light km/s
dMat = rEst - rSats;
dMag = vecnorm(dMat);
dMatNorm = dMat./dMag;
zEst = 1/c*(dMag(2:end) - dMag(1)).';
hMatEst = 1/c*(dMatNorm(:,2:end) - dMatNorm(:,1)).';
end