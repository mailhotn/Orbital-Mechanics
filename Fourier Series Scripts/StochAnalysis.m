clear
%%
dataFolder = 'C:\Users\User\Google Drive\Doc Data\Stochastic Estimation';
% load([dataFolder '\name.mat']);
% load([dataFolder '\StochErr_12-7-2023_15-10.mat']);
% load([dataFolder '\StochErr_12-7-2023_17-4.mat']);
% load([dataFolder '\StochErr_16-7-2023_14-41.mat']);
load([dataFolder '\StochErr_19-7-2023_17-0.mat']); % 20,000 runs, no wrapping

%% Analyze Data
t = errData.t;
nT = length(t);
sigB = nan(6,nT);
sigF = nan(6,nT);
for iOe = 1:6
    sigB(iOe,:) = std(squeeze(errData.errorB(iOe,:,:)));
    sigF(iOe,:) = std(squeeze(errData.errorF(iOe,:,:)));
end
mErrB = squeeze(mean(errData.errorB,2));
mErrF = squeeze(mean(errData.errorF,2));

%% Plot STD

figure(1)
plot(t,sigB(1,:),'o',t,sigF(1,:),'.')

figure(2)
plot(t,sigB(2,:),'o',t,sigF(2,:),'.')

figure(3)
plot(t,sigB(3,:),'o',t,sigF(3,:),'.')

figure(4)
plot(t,sigB(4,:),'o',t,sigF(4,:),'.')

figure(5)
plot(t,sigB(5,:),'o',t,sigF(5,:),'.')

figure(6)
plot(t,sigB(6,:),'o',t,sigF(6,:),'.')

%% Plot Mean Errors
figure(11)
plot(t,mErrB(1,:),'o',t,mErrF(1,:),'.')

figure(12)
plot(t,mErrB(2,:),'o',t,mErrF(2,:),'.')

figure(13)
plot(t,mErrB(3,:),'o',t,mErrF(3,:),'.')

figure(14)
plot(t,mErrB(4,:),'o',t,mErrF(4,:),'.')

figure(15)
plot(t,mErrB(5,:),'o',t,mErrF(5,:),'.')

figure(16)
plot(t,mErrB(6,:),'o',t,mErrF(6,:),'.')
