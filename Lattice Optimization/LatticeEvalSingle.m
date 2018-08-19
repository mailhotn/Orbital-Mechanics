
datafolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\Lattice Version 1';
walkerFolder = ['C:\Users\User\Dropbox\Walker Optimization Data'...
    '\Previous Optimization Runs\Walker RGT Ex Search delta inc 10'];
%% Choose Constellation
latGs = 50;
lonGs = 0;
ecc = 0.05;
nSats = 51;
nPlanes = 51;
nC2 = 1;
Arch.nDays = 1;
Arch.nRepeats = 14; 
load ([datafolder '\OptParams.mat']);
load([datafolder '\LatticeExSol_Lat_' num2str(latGs) ...
    '_nSats_' num2str(nSats) '_ecc_' num2str(ecc) '.mat']);
%% Create Constellation
Arch.nSats = nSats;
Arch.nPlanes = nPlanes;
Phase.nC3 = Arch.nPlanes; % Constraint! can be different, would increase number of relative orbits.
Arch.nSatsPerAop = Arch.nDays; % Constraint!  can be different if nDays > 1
Arch.nAops = Arch.nSats/Arch.nPlanes/Arch.nSatsPerAop;
gcdOrbits = gcd(Arch.nPlanes,Phase.nC3);
if gcdOrbits > Arch.nRepeats
    Phase.nC1 = Arch.nRepeats;
else
    l = floor(Arch.nRepeats/gcdOrbits);
    if l*gcdOrbits == Arch.nRepeats
        l = l-1;
    end
    Phase.nC1 = Arch.nRepeats - l*gcdOrbits;
end
Phase.nC2 = nC2;

LC = LatticeConstellation(Arch,Phase,ExSol.Orbit,ExSol.InitCon);
nSatsW = 64;
load([walkerFolder '\WalkerRgtExSol_Lat_' num2str(latGs)...
                '_T_' num2str(nSatsW) '.mat']);
WC = WalkerConstellation(64,8,2,ExSol.inc,ExSol.alt,ExSol.raan0);
%% Propagate
Prop = Propagator(LC,PropParams.relTol,PropParams.absTol);
PropW = Propagator(WC,PropParams.relTol,PropParams.absTol);
[propTime, propState] = Prop.PropEciJ2(PropParams.timeVec);
[propTimeW, propStateW] = PropW.PropEciJ2(PropParams.timeVec);
% Evaluate PDOP
[pdop, satsIs] = TdoaPdopVec(propState,propTime,latGs,0,0,PropParams.elevMin);
%% Plot Stuff
figure(1)
PlotGroundTrack(propState,propTime,0)

% figure(2)
% yyaxis left
% plot(propTime/3600,pdop)
% ylabel('PDOP')
% yyaxis right
% plot(propTime/3600,satsIs)
% ylabel('Sats In Sight')
% xlabel('Time [hr]')
% grid
%% Plot PDOP Map
lats = 40:1:60;
lons = -20:1:20;
[LON,LAT] = meshgrid(lons,lats);
PDOP = nan(size(LON));
PDOPW = nan(size(LON));
tic
parfor iLon = 1:length(lons)
    PDOP2 = zeros(length(lats),1);
    PDOPW2 = zeros(length(lats),1);
    for iLat = 1:length(lats)
        pdop = TdoaPdopVec(propState,propTime,LAT(iLat,iLon),LON(iLat,iLon)...
            ,0,PropParams.elevMin);
        pdop(pdop>100) = 100;
        pdop(isnan(pdop)) = 100;
        PDOP2(iLat) = mean(pdop(~isnan(pdop)));
        
        pdopW = TdoaPdopVec(propStateW,propTimeW,LAT(iLat,iLon),LON(iLat,iLon)...
            ,0,PropParams.elevMin);
        pdopW(pdopW>100) = 100;
        pdopW(isnan(pdopW)) = 100;
        PDOPW2(iLat) = mean(pdopW(~isnan(pdopW)));
    end
    PDOP(:,iLon) =PDOP2;
    PDOPW(:,iLon) =PDOPW2;
end
toc
figure(3)
subplot(2,1,1)

contourf(LON,LAT,PDOP,2000,'LineColor','none')
title('Lattice Flower Constellation 51/51')
hold on

colormap jet
shading interp
c = colorbar;
c.Label.String = 'Mean PDOP';
xlabel('Longitude')
ylabel('Latitude')
load coastlines
geoshow(coastlat,coastlon)
plot(lonGs,latGs,'y*','LineWidth',1.5)
axis equal
xlim([-30,30])
ylim([30,70])
grid on
hold off

subplot(2,1,2)

contourf(LON,LAT,PDOPW,2000,'LineColor','none')
title('Walker Constellation 64/8')
hold on

colormap jet
shading interp
c = colorbar;
c.Label.String = 'Mean PDOP';
xlabel('Longitude')
ylabel('Latitude')
load coastlines
geoshow(coastlat,coastlon)
plot(lonGs,latGs,'y*','LineWidth',1.5)
axis equal
xlim([-30,30])
ylim([30,70])
grid on
hold off