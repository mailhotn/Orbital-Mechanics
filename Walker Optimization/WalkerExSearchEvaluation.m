datafolder = 'C:\Users\User\Dropbox\Walker Optimization Data';

maxSats = 100;
minSats = 20;
nCons = maxSats - minSats + 1;
latList = 10:10:80;
nLats = length(latList);
maxPdop   = nan(nLats,nCons);
maxPdop2  = nan(nLats,nCons);
intPdop   = nan(nLats,nCons);
intPdop2  = nan(nLats,nCons);
meanPdop  = nan(nLats,nCons);
coverage  = nan(nLats,nCons);
inc       = nan(nLats,nCons);
alt       = nan(nLats,nCons);
nSatsT    = repmat(minSats:maxSats,nLats,1);

%%
for iLat = 1:length(latList)
    for iSat = 1:nCons
        load([datafolder '\WalkerRgtExSol_Lat_' num2str(latList(iLat))...
            '_T_' num2str(nSatsT(iLat,iSat)) '.mat']);
        maxPdop(iLat,iSat) = ExSol.maxPdop(ExSol.optF+1,ExSol.optP);
        intPdop(iLat,iSat) = ExSol.intPdop(ExSol.optF+1,ExSol.optP);
        meanPdop(iLat,iSat) = ExSol.meanPdop(ExSol.optF+1,ExSol.optP);
        coverage(iLat,iSat) = ExSol.coverage(ExSol.optF+1,ExSol.optP);
        inc(iLat,iSat) = ExSol.inc;
        alt(iLat,iSat) = ExSol.alt;
        % Find Optimal Point for Max PDOP
        
        [minForP,indP] = min(sqrt(ExSol.maxPdop.^2 + ExSol.intPdop.^2));
        [~,optP] = min(minForP);
        optF = indP(optP) - 1;
        maxPdop2(iLat,iSat) = ExSol.maxPdop(optF+1,optP);
        intPdop2(iLat,iSat) = ExSol.intPdop(optF+1,optP);
    end
end

figure(1) % integral PDOP
semilogy(nSatsT(1:7,:).',intPdop(1:7,:).','o',...
    nSatsT(8,:),intPdop(8,:),'s')
title('Integral of PDOP')
legend('\phi_0=10','\phi_0=20','\phi_0=30','\phi_0=40','\phi_0=50','\phi_0=60',...
    '\phi_0=70','\phi_0=80')
xlabel('# Sats')
ylabel('$\frac{1}{T}\int^{T}_{0}{PDOPdt}$','interpreter','latex','fontsize',12)
% ylim([0,10])
grid minor

figure(2) % max PDOP
semilogy(nSatsT(1:7,:).',maxPdop(1:7,:).','o',...
    nSatsT(8,:),maxPdop(8,:),'s')
title('Max PDOP')
legend('\phi_0=10','\phi_0=20','\phi_0=30','\phi_0=40','\phi_0=50','\phi_0=60',...
    '\phi_0=70','\phi_0=80')
xlabel('# Sats')
ylabel('$\max{PDOP}$','interpreter','latex')
% ylim([0,10000])
grid minor


figure(3) % max PDOP as cost
semilogy(nSatsT(1:7,:).',maxPdop2(1:7,:).','o',...
    nSatsT(8,:),maxPdop2(8,:),'s')
title('Max PDOP')
legend('\phi_0=10','\phi_0=20','\phi_0=30','\phi_0=40','\phi_0=50','\phi_0=60',...
    '\phi_0=70','\phi_0=80')
xlabel('# Sats')
ylabel('$\max{PDOP}$','interpreter','latex')
% ylim([0,10000])
grid minor


figure(4) % integral PDOP
semilogy(nSatsT(1:7,:).',intPdop2(1:7,:).','o',...
    nSatsT(8,:),intPdop2(8,:),'s')
title('Integral of PDOP')
legend('\phi_0=10','\phi_0=20','\phi_0=30','\phi_0=40','\phi_0=50','\phi_0=60',...
    '\phi_0=70','\phi_0=80')
xlabel('# Sats')
ylabel('$\frac{1}{T}\int^{T}_{0}{PDOPdt}$','interpreter','latex','fontsize',12)
% ylim([0,10])
grid minor

figure(5) % coverage
plot(nSatsT(1:7,:).',coverage(1:7,:).','o',...
    nSatsT(8,:),coverage(8,:),'s')
title('Coverage')
ylabel('Coverage %')
legend('\phi_0=10','\phi_0=20','\phi_0=30','\phi_0=40','\phi_0=50','\phi_0=60',...
    '\phi_0=70','\phi_0=80')
xlabel('# Sats')
grid minor
ylim([0,100])

figure(6) % Altitude
plot(nSatsT(1:7,:).',alt(1:7,:).','o',...
        nSatsT(8,:),alt(8,:),'s')
title('Altitude')
legend('\phi_0=10','\phi_0=20','\phi_0=30','\phi_0=40','\phi_0=50','\phi_0=60',...
    '\phi_0=70','\phi_0=80')
% ylim([0,100])
xlabel('# Sats')
ylabel('Altitude [km]')
grid minor

figure(7) % Inclination
plot(nSatsT(1:7,:).',inc(1:7,:).','o',...
        nSatsT(8,:),inc(8,:),'s')
title('Inclination')
legend('\phi_0=10','\phi_0=20','\phi_0=30','\phi_0=40','\phi_0=50','\phi_0=60',...
    '\phi_0=70','\phi_0=80')
% ylim([0,100])
xlabel('# Sats')
ylabel('Inclination [°]')
grid minor