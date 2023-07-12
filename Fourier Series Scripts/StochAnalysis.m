dataFolder = 'C:\Users\User\Dropbox\Doc Fourier Data\Stochastic Estimation';
% load([dataFolder '\name.mat']);

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
