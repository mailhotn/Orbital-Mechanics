%% Define Constellation & Propagate
%
c = 3e5;
timeVar = (1e-6)^2;
datafolder = 'C:\Users\User\Google Drive\Master''s Degree\Lattice Optimization Data';
latGs = 50;
lonGs = 0;
nSats = 80;
ecc = 0.02;
load ([datafolder '\OptParams.mat']);
load([datafolder '\LatticeExSol_Lat_' num2str(latGs) ...
    '_nSats_' num2str(nSats) '_ecc_' num2str(ecc) '.mat']);
Arch.nDays = 1;
Arch.nRepeats = 14;
Arch.nSats = nSats;
Arch.nPlanes = ExSol.optNPlanes;
iPlane = ExSol.archMat(1,:) == Arch.nPlanes;
Arch.nAops = ExSol.archMat(2,iPlane);
Arch.nSatsPerAop = ExSol.archMat(3,iPlane);

Phase.nC1 = ExSol.phaseMat(1,iPlane);
Phase.nC2 = ExSol.phaseMat(2,iPlane);
Phase.nC3 = ExSol.phaseMat(3,iPlane);

LC = LatticeConstellation(Arch,Phase,ExSol.Orbit,ExSol.InitCon);

Prop = Propagator(LC,PropParams.relTol,PropParams.absTol);
[time, xEci] = Prop.PropEciJ2(PropParams.timeVec);
%% 3 Satellite Spherical Earth
tic
[rGsEst, estError, pdop, flagVec] = Geolocate3SatsVec(xEci, time, latGs, lonGs,...
    PropParams.elevMin, timeVar);
toc
% Plot PDOP & Errors
figure(1)
semilogy(time/3600,estError,'.',time/3600,pdop*c*sqrt(timeVar))
xlim([0,24])
xlabel('Time [hr]')
legend('Geolocation Error [km]','Ambiguous Solution Error [km]'...
    ,'PDOP\cdotc\cdot\sigma_t','location','best')
% Plot PDOP & Flags
isGoodSol  = flagVec <= 3;
pdopGood = pdop;
pdopGood(~isGoodSol) = nan;

isNoSol    = flagVec == 4;
pdopNo  = pdop;
pdopNo(~isNoSol) = nan;

isAmbigSol = flagVec == 5;
pdopAmbig = pdop;
pdopAmbig(~isAmbigSol) = nan;

figure(2)
semilogy(time/3600,pdopGood,'.',time/3600,pdopNo,'.',...
    time/3600,pdopAmbig,'.')
xlim([0 24])
legend('Unambiguous Solution','No Solution','Ambiguous Solution')
xlabel('Time [hr]')
ylabel('PDOP')

disp([newline newline 'Time Variance: (' num2str(sqrt(timeVar)) ')^2 sec^2' newline...
      'Unambiguous Ratio: ' num2str(sum(flagVec<=3)/numel(flagVec)*100,4) '%' newline...
      'Ambiguous Ratio: ' num2str(sum(flagVec==5)/numel(flagVec)*100,4) '%' newline...
      'No Solution Ratio: ' num2str(sum(flagVec==4)/numel(flagVec)*100,4) '%' newline...
      'Mean Unambiguous PDOP: ' num2str(mean(pdopGood(~isnan(pdopGood))),4)...
      newline 'Mean Ambiguous PDOP: ' num2str(mean(pdopAmbig(~isnan(pdopAmbig))),4)...
      newline 'Mean No Solution PDOP: ' num2str(mean(pdopNo(~isnan(pdopNo))),4)])
%% N Satellite Spherical Earth
tic
[rGsEst, estError, pdop, flagVec] = GeolocateSphereVec(xEci, time, latGs, lonGs,...
    PropParams.elevMin, timeVar);
toc
% Plot PDOP & Errors
figure(1)
semilogy(time/3600,estError,'.')
xlim([0,24])
xlabel('Time [hr]')
legend('Geolocation Error [km]','Ambiguous Solution Error [km]'...
    ,'PDOP\cdotc\cdot\sigma_t','location','best')