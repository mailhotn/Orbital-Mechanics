datafolder = 'C:\Users\User\Dropbox\Walker Optimization Data\Previous Optimization Runs\';
% load([datafolder '\OptParams.mat']);
maxSats = 80;
minSats = 20;
nCons = maxSats - minSats + 1;
latList = 10:10:80;
nLats = length(latList);
nSatsT    = repmat(minSats:maxSats,nLats,1);
incList = [0,5,10,15,20,25];

intTarget = 2;
maxTarget = 10;
nSatsToAchieve = nan(nLats,length(incList));
nPlanesToAchieve = nan(nLats,length(incList));

%%
for iInc = 1:length(incList)
    for iLat = 1:length(latList)
        for iSat = 1:nCons
            load([datafolder 'Walker RGT Ex Search delta inc ' ...
                num2str(incList(iInc)) '\WalkerRgtExSol_Lat_' ...
                num2str(latList(iLat)) '_T_' num2str(nSatsT(iLat,iSat)) '.mat']);
            if isnan(nSatsToAchieve(iLat,iInc))
                achieveInt = ExSol.intPdop < intTarget;
                achieveMax = ExSol.maxPdop < maxTarget;
                if any(any(achieveInt & achieveMax))
                    nSatsToAchieve(iLat,iInc) = nSatsT(iLat,iSat);
                    planeVec = sum(achieveInt & achieveMax,1);
                    [~,nPlanesToAchieve(iLat,iInc)] = min(~(planeVec>0)); 
                end
            end
        end
    end
end
[minT,indMin] = min(nSatsToAchieve.');
figure(1)
plot(latList,incList(indMin),'o')
xlabel('Ground Station Latitude [°]')
ylabel('\Deltai [°]')

figure(2)
plot(incList,nSatsToAchieve(1:7,:),...
     incList,nSatsToAchieve(8,:),'--','Linewidth',2)
xlabel('\Deltai [°]')
ylabel('# Satellites')
legend('\phi_0=10','\phi_0=20','\phi_0=30','\phi_0=40','\phi_0=50','\phi_0=60',...
    '\phi_0=70','\phi_0=80')