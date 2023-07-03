clear
%% Define Satellite & Noise Model
ic = [7000, 0.01, 55, 10, 10, 10];
primary = earth();
mu = primary.mu;
Sat = SingleSat(ic);
Prop = Propagator(Sat);

sigR = 5e-3;
sigV = 2e-5;
covEci = diag([sigR*ones(1,3),sigV*ones(1,3)]);
kMax = 5;

t = 0:100:86400;
nT = length(t);

nMonte = 10000;

errorB = nan(6,nMonte,length(t));
errorF = nan(6,nMonte,length(t));
%% True Prop

[~, eciTrue] = Prop.PropEciJ2(t);
oeTrue = eci2oe(eciTrue);
mOeTrue = osc2me(oeTrue);
mOeTrue(5:6,:) = wrapTo180(mOeTrue(5:6,:));

%% Monte Carlo Runs
parfor iMonte = 1:nMonte
    %% Contaminate
noise = mvnrnd(zeros(6,1),covEci,nT);
eciMeas = eciTrue + noise;
oeMeas = eci2oe(eciMeas);

%% Calculate Means

mOeMeas = osc2me(oeMeas);
mOeMeas(5:6,:) = wrapTo180(mOeMeas(5:6,:));
mOeFour = nan(size(oeTrue));

[~, lpeSpec] = LpeJ2Fourier(oeMeas,kMax,primary);
M = oeMeas(6,:);
a = oeMeas(1,:);
nMo = sqrt(mu./a.^3);
k = (1:kMax).';
Ck = -cos(k*M)./k./nMo;
Sk = sin(k*M)./k./nMo;

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

mOeFour = mOeMeas - [smaVar;eccVar;incVar;ranVar;aopVar;manVar];
% for iTime = 1:nT
%     [~, lpeSpec] = LpeJ2Fourier(oeMeas(:,iTime),kMax,primary);
%     k = 1:kMax;
%     M = oeMeas(6,iTime);
%     a = oeMeas(1,iTime);
%     nMo = sqrt(mu/a^3);
%     % find M
%     trigVec1 = [cos(k*M)./k./nMo,sin(k*M)./k/nMo].';
%     M2 = oeMeas(6,iTime) - lpeSpec(6,:)*trigVec1;
%     % Use updated M for estimate
%     trigVec2 = [cos(k*M2)./k/nMo,sin(k*M2)./k/nMo].';
%     mOeFour(:,iTime) = oeMeas(:,iTime) - lpeSpec*trigVec2;
% end
mOeFour(5:6,:) = wrapTo180(mOeFour(5:6,:));
errorB(:,iMonte,:) = mOeMeas - mOeTrue;
errorF(:,iMonte,:) = mOeFour - mOeTrue;
end
sigB = nan(6,nT);
sigF = nan(6,nT);
for iOe = 1:6
    sigB(iOe,:) = std(squeeze(errorB(iOe,:,:)));
    sigF(iOe,:) = std(squeeze(errorF(iOe,:,:)));
end

%% Plot STD

figure(1)
plot(t,sigB(1,:),'-',t,sigF(1,:),'.')

figure(2)
plot(t,sigB(2,:),'-',t,sigF(2,:),'.')

figure(3)
plot(t,sigB(3,:),'-',t,sigF(3,:),'.')

figure(4)
plot(t,sigB(4,:),'-',t,sigF(4,:),'.')

figure(5)
plot(t,sigB(5,:),'-',t,sigF(5,:),'.')

figure(6)
plot(t,sigB(6,:),'-',t,sigF(6,:),'.')


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
