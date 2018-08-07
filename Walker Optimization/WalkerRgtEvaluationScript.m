%% Initializaton
% Simulation Parameters
timeVector = 0:10:86164;
lonGs = 0;
elevMin  = 5;
relTol = 1e-6;
absTol = 1e-6;

% Initialize data vectors
minSats = 20;
maxSats  = 80;
nCons = maxSats - minSats + 1;
minLat = 10;
maxLat = 80;
nLats = length(minLat:10:maxLat);
primeList = primes(maxSats);
latList   = minLat:10:maxLat;
nSatsT    = repmat(minSats:maxSats,nLats,1);
nPlanesP  = nan(nLats,nCons);
phasingF  = nan(nLats,nCons);
gaInc     = nan(nLats,nCons);
gaAlt     = nan(nLats,nCons);
gaGmst0   = nan(nLats,nCons);
maxPdop   = nan(nLats,nCons);
intPdop   = nan(nLats,nCons);
meanPdop  = nan(nLats,nCons);
coverage  = nan(nLats,nCons);

datafolder = 'C:\Users\User\Dropbox\Walker Optimization Data\Previous Optimization Runs\Walker RGT GA mean pdop';

%% Go over all solutions & sim Osculating to get new fitness
tic
for iLat = 1:nLats
    for iSat = 1:nCons
        if ~any(nSatsT(iLat,iSat) == primeList)
            % Load Data & Assign
            load([datafolder '\WalkerMeanRgtSolLat_' num2str(latList(iLat)) ...
                '_T_' num2str(nSatsT(iLat,iSat)) '.mat']);
            nPlanesP(iLat,iSat) = GaRgtSol.nPlanesP;
            phasingF(iLat,iSat) = GaRgtSol.phasingF;
            gaInc(iLat,iSat)    = GaRgtSol.inc;
            gaAlt(iLat,iSat)    = CalcRgtElement([],0,GaRgtSol.inc,GaRgtSol.jRepeats,1)-6378.137;
            % Simulate Constellation
            WC = WalkerConstellation(GaRgtSol.nSatsT,GaRgtSol.nPlanesP,GaRgtSol.phasingF,...
                GaRgtSol.inc,gaAlt(iLat,iSat));
            Prop = Propagator(WC,relTol,absTol);
            [propTime, propState] = Prop.PropEciJ2(timeVector);
            [pdop, satsIs] = TdoaPdopVec(propState,propTime,GaRgtSol.latGs,GaRgtSol.lonGs,...
                GaRgtSol.gmst0, GaRgtSol.elevMin);
            % Calculate Performance Indeces
            coverage(iLat,iSat) = 100 - sum(isnan(pdop))/length(pdop)*100;
            maxPdop(iLat,iSat) = max(pdop);
            pdop(pdop>1000) = 1000;
            pdop(isnan(pdop)) = 1000;
            meanPdop(iLat,iSat) = mean(pdop(~isnan(pdop)));
%             pdop(isnan(pdop)) = max(pdop)*10;
            intPdop(iLat,iSat) = trapz(propTime,pdop)/(propTime(end)-propTime(1));
        end
    end
end
toc
%% Plot Results
figure(1) % integral PDOP
semilogy(nSatsT(1,:),intPdop(1,:),'o',...
    nSatsT(2,:),intPdop(2,:),'o',...
    nSatsT(3,:),intPdop(3,:),'o',...
    nSatsT(4,:),intPdop(4,:),'o',...
    nSatsT(5,:),intPdop(5,:),'o',...
    nSatsT(6,:),intPdop(6,:),'o',...
    nSatsT(7,:),intPdop(7,:),'o',...
    nSatsT(8,:),intPdop(8,:),'s')
title('Integral of PDOP')
legend('\phi_0=10','\phi_0=20','\phi_0=30','\phi_0=40','\phi_0=50','\phi_0=60',...
    '\phi_0=70','\phi_0=80')
xlabel('# Sats')
ylabel('$\frac{1}{T}\int^{T}_{0}{PDOPdt}$','interpreter','latex','fontsize',12)
% ylim([0,10])
grid minor

figure(2) % max PDOP
semilogy(nSatsT(1,:),maxPdop(1,:),'o',...
    nSatsT(2,:),maxPdop(2,:),'o',...
    nSatsT(3,:),maxPdop(3,:),'o',...
    nSatsT(4,:),maxPdop(4,:),'o',...
    nSatsT(5,:),maxPdop(5,:),'o',...
    nSatsT(6,:),maxPdop(6,:),'o',...
    nSatsT(7,:),maxPdop(7,:),'o',...
    nSatsT(8,:),maxPdop(8,:),'s')
title('Max PDOP')
legend('\phi_0=10','\phi_0=20','\phi_0=30','\phi_0=40','\phi_0=50','\phi_0=60',...
    '\phi_0=70','\phi_0=80')
xlabel('# Sats')
ylabel('$\max{PDOP}$','interpreter','latex')
% ylim([0,10000])
grid minor

figure(3) % mean PDOP
semilogy(nSatsT(1,:),meanPdop(1,:),'o',...
    nSatsT(2,:),meanPdop(2,:),'o',...
    nSatsT(3,:),meanPdop(3,:),'o',...
    nSatsT(4,:),meanPdop(4,:),'o',...
    nSatsT(5,:),meanPdop(5,:),'o',...
    nSatsT(6,:),meanPdop(6,:),'o',...
    nSatsT(7,:),meanPdop(7,:),'o',...
    nSatsT(8,:),meanPdop(8,:),'s')
title('Mean PDOP')
legend('\phi_0=10','\phi_0=20','\phi_0=30','\phi_0=40','\phi_0=50','\phi_0=60',...
    '\phi_0=70','\phi_0=80')
xlabel('# Sats')
ylabel('mean PDOP')
grid minor
% ylim([0,100])

figure(4) % coverage
plot(nSatsT(1,:),coverage(1,:),'o',...
    nSatsT(2,:),coverage(2,:),'o',...
    nSatsT(3,:),coverage(3,:),'o',...
    nSatsT(4,:),coverage(4,:),'o',...
    nSatsT(5,:),coverage(5,:),'o',...
    nSatsT(6,:),coverage(6,:),'o',...
    nSatsT(7,:),coverage(7,:),'o',...
    nSatsT(8,:),coverage(8,:),'s')
title('Coverage')
ylabel('Coverage %')
legend('\phi_0=10','\phi_0=20','\phi_0=30','\phi_0=40','\phi_0=50','\phi_0=60',...
    '\phi_0=70','\phi_0=80')
xlabel('# Sats')
grid minor
ylim([0,100])

figure(5) % Altitude
plot(nSatsT(1,:),gaAlt(1,:),'o',...
    nSatsT(2,:),gaAlt(2,:),'o',...
    nSatsT(3,:),gaAlt(3,:),'o',...
    nSatsT(4,:),gaAlt(4,:),'o',...
    nSatsT(5,:),gaAlt(5,:),'o',...
    nSatsT(6,:),gaAlt(6,:),'o',...
    nSatsT(7,:),gaAlt(7,:),'o',...
    nSatsT(8,:),gaAlt(8,:),'s')
title('Altitude')
legend('\phi_0=10','\phi_0=20','\phi_0=30','\phi_0=40','\phi_0=50','\phi_0=60',...
    '\phi_0=70','\phi_0=80')
% ylim([0,100])
xlabel('# Sats')
ylabel('Altitude [km]')
grid minor

figure(6) % Inclination
plot(nSatsT(1:7,:).',gaInc(1:7,:),'o',...
    nSatsT(8,:),gaInc(8,:),'s')
title('Inclination')
legend('\phi_0=10','\phi_0=20','\phi_0=30','\phi_0=40','\phi_0=50','\phi_0=60',...
    '\phi_0=70','\phi_0=80')
% ylim([0,100])
xlabel('# Sats')
ylabel('Inclination [°]')
grid minor
%% Save Data
GaEval.nSatsT   = nSatsT;
GaEval.nPlanesP = nPlanesP;
GaEval.phasingF = phasingF;
GaEval.inc      = gaInc;
GaEval.alt      = gaAlt;
GaEval.latList  = latList;
GaEval.gmst0    = gaGmst0;
GaEval.maxPdop  = maxPdop;
GaEval.intPdop  = intPdop;
GaEval.meanPdop = meanPdop;
GaEval.coverage = coverage;
GaEval.elevMin  = elevMin;
GaEval.lonGs    = lonGs;
GaEval.relTol   = relTol;
GaEval.absTol   = absTol;

save([datafolder '\GA Evaluation Data.mat'],'GaEval');
