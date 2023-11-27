clear
%%
dataFolder = 'C:\Users\User\Google Drive\Doc Data\Stochastic Estimation';
% load([dataFolder '\name.mat']);
% load([dataFolder '\StochErr_12-7-2023_15-10.mat']);
% load([dataFolder '\StochErr_12-7-2023_17-4.mat']);
% load([dataFolder '\StochErr_16-7-2023_14-41.mat']);
% load([dataFolder '\StochErr_19-7-2023_17-0.mat']); % 20,000 runs, no wrapping

% Added element saving
% load([dataFolder '\StochErr_7-8-2023_14-31.mat']); % 20,000 runs, no wrapping

% Nonrandom Orbit
% load([dataFolder '\StochErr_9-8-2023_12-53.mat']); 

% Back to Random, use nominal osc sma for mean motion in Fourier
% load([dataFolder '\StochErr_16-8-2023_15-49.mat']); 

% nominal mean sma for mean motion
load([dataFolder '\StochErr_17-8-2023_10-45.mat']);
% load([dataFolder '\StochErr_17-8-2023_11-15.mat']); % long period

% Fixed(?) error with Fourier using mean? - or not
% load([dataFolder '\StochErr_27-11-2023_13-32.mat']); 





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

%% Stats
meanSigF = mean(sigF,2)
meanSigB = mean(sigB,2)
sigDiff = meanSigB - meanSigF
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
