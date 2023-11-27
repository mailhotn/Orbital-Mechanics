clear
%% Define Satellite & Noise Model
dataFolder = 'C:\Users\User\Google Drive\Doc Data\Stochastic Estimation';
drivePath = 'C:\Users\User\Google Drive'; % ASRI
smaRange = [7000,8000];
eccRange = [-2,-1]; % 10^x
incRange = [1,30];
primary = earth();
mu = primary.mu;

sigR = 5e-3;
sigV = 2e-5;
covEci = diag([sigR*ones(1,3),sigV*ones(1,3)]);
kMax = 5;

t = 0:100:86400;
nT = length(t);

nMonte = 1000;

errorB = nan(6,nMonte,length(t));
errorF = nan(6,nMonte,length(t));
oeErr = nan(6,nMonte,length(t));
mOeMat = nan(6,nMonte,length(t));
%% Monte Carlo Runs
totalTime = tic;
pool = parpool('Threads');
parfor iMonte = 1:nMonte
    % General IC
    % sma = smaRange(1) + rand*(smaRange(2)-smaRange(1));
    % inc = incRange(1) + rand*(incRange(2)-incRange(1));
    % ecc = 10^(eccRange(1) + rand*(eccRange(2)-eccRange(1)));
    % ic = [sma, ecc, inc, rand(1,3)*360].';
    % Fixed IC
    sma = 7000;
    ecc = 0.01;
    inc = 30;
    ic = [sma, ecc, inc, 10, 10, 10].';
    % Initialize Constellation
    icM = osc2me(ic);
    Sat = SingleSat(ic);
    Prop = Propagator(Sat);
    %% True Prop
    [~, eciTrue] = Prop.PropEciJ2(t);
    oeTrue = eci2oe(eciTrue);
    mOeTrue = osc2me(oeTrue);
    mOeMat(:,iMonte,:) = mOeTrue;
    % mOeTrue(5:6,:) = wrapTo180(mOeTrue(5:6,:));
    %% Contaminate
    noise = mvnrnd(zeros(6,1),covEci,nT);
    eciMeas = eciTrue + noise;
    oeMeas = eci2oe(eciMeas);
    oeErr(:,iMonte,:) = oeMeas - oeTrue;
    %% Calculate Means
    mOeMeas = osc2me(oeMeas);
    % mOeMeas(4:6,:) = wrapTo180(mOeMeas(4:6,:));
    mOeFour = nan(size(oeTrue));

    [~, lpeSpec] = LpeJ2Fourier(oeMeas,kMax,primary);
    M = oeMeas(6,:);
    % a = oeMeas(1,:); % Use nominal a instead - cheating? not really
    nomSma = icM(1);
    nMo = sqrt(mu./nomSma.^3);
    k = (1:kMax).';
    Ck = -cos(k*M)./k./nMo;
    Sk = sin(k*M)./k./nMo;

    % Fix M
    manVar =  sum(Sk.*squeeze(lpeSpec(11,:,:)).',1) + ...
        sum(Ck.*squeeze(lpeSpec(12,:,:)).',1);
    M2 = oeMeas(6,:)- manVar;
    Ck = -cos(k*M2)./k./nMo;
    Sk = sin(k*M2)./k./nMo;
    % Calculate other Elements
    smaVar =  sum(Sk.*squeeze(lpeSpec(1,:,:)).',1) + ...
        sum(Ck.*squeeze(lpeSpec(2,:,:)).',1);
    eccVar =  sum(Sk.*squeeze(lpeSpec(3,:,:)).',1) + ...
        sum(Ck.*squeeze(lpeSpec(4,:,:)).',1);
    incVar =  sum(Sk.*squeeze(lpeSpec(5,:,:)).',1) + ...
        sum(Ck.*squeeze(lpeSpec(6,:,:)).',1);
    ranVar =  sum(Sk.*squeeze(lpeSpec(7,:,:)).',1) + ...
        sum(Ck.*squeeze(lpeSpec(8,:,:)).',1);
    aopVar =  sum(Sk.*squeeze(lpeSpec(9,:,:)).',1) + ...
        sum(Ck.*squeeze(lpeSpec(10,:,:)).',1);
    manVar =  sum(Sk.*squeeze(lpeSpec(11,:,:)).',1) + ...
        sum(Ck.*squeeze(lpeSpec(12,:,:)).',1);

    mOeFour = oeMeas - [smaVar;eccVar;incVar;ranVar;aopVar;manVar];
    % mOeTrue(4:6,:) = wrapTo180(mOeTrue(4:6,:));
    errorB(:,iMonte,:) = mOeMeas - mOeTrue;
    errorF(:,iMonte,:) = mOeFour - mOeTrue;
end
delete(pool);
sigOe = nan(6,nT);
for iOe = 1:6
    sigOe(iOe,:) = std(squeeze(oeErr(iOe,:,:)));
end
covOe = mean(sigOe,2);
eTime = toc(totalTime);
reportIFTTT(drivePath,eTime);
%% Save Data
errData = struct();
errData.smaRange = smaRange;
errData.eccRange = eccRange;
errData.incRange = incRange;
errData.kMax = kMax;
errData.nMonte = nMonte;
errData.covEci = covEci;
errData.t = t;
errData.runTime = eTime;
errData.errorB = single(errorB);
errData.errorF = single(errorF);
errData.mOeTrue = single(mOeMat);
errDat.covOe = covOe;

c = clock;
save([dataFolder '\StochErr_' num2str(c(3)) '-' num2str(c(2)) '-' ...
    num2str(c(1)) '_' num2str(c(4)) '-' num2str(c(5)), '.mat'],'errData');

% sigB = nan(6,nT);
% sigF = nan(6,nT);
% for iOe = 1:6
%     sigB(iOe,:) = std(squeeze(errorB(iOe,:,:)));
%     sigF(iOe,:) = std(squeeze(errorF(iOe,:,:)));
% end
% 
% %% Plot STD
% 
% figure(1)
% plot(t,sigB(1,:),'-',t,sigF(1,:),'.')
% 
% figure(2)
% plot(t,sigB(2,:),'-',t,sigF(2,:),'.')
% 
% figure(3)
% plot(t,sigB(3,:),'-',t,sigF(3,:),'.')
% 
% figure(4)
% plot(t,sigB(4,:),'-',t,sigF(4,:),'.')
% 
% figure(5)
% plot(t,sigB(5,:),'-',t,sigF(5,:),'.')
% 
% figure(6)
% plot(t,sigB(6,:),'-',t,sigF(6,:),'.')


% %% Plot Single Run
%
% figure(1)
% plot(t,mOeMeas(1,:),t,mOeTrue(1,:),t,mOeFour(1,:))
%
% figure(2)
% plot(t,mOeMeas(2,:),t,mOeTrue(2,:),t,mOeFour(2,:))
%
% figure(3)
% plot(t,mOeMeas(3,:),t,mOeTrue(3,:),t,mOeFour(3,:))
%
% figure(4)
% plot(t,mOeMeas(4,:),t,mOeTrue(4,:),t,mOeFour(4,:))
%
% figure(5)
% plot(t,mOeMeas(5,:),t,mOeTrue(5,:),t,mOeFour(5,:))
%
% figure(6)
% plot(t,mOeMeas(6,:),t,mOeTrue(6,:),t,mOeFour(6,:))
%
