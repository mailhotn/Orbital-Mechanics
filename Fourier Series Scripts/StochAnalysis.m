dataFolder = 'C:\Users\User\Google Drive\Doc Data\Stochastic Estimation';
% load([dataFolder '\name.mat']);
% load([dataFolder '\StochErr_12-7-2023_15-10.mat']);
% load([dataFolder '\StochErr_12-7-2023_17-4.mat']);
% load([dataFolder '\StochErr_16-7-2023_14-41.mat']);
load([dataFolder '\StochErr_30-7-2023_14-45.mat']);

%% Analyze Data
t = errData.t;
nT = length(t);
sigB = nan(6,nT);
sigF = nan(6,nT);
for iOe = 1:6
    sigB(iOe,:) = std(squeeze(errData.errorB(iOe,:,:)));
    sigF(iOe,:) = std(squeeze(errData.errorF(iOe,:,:)));
end

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
