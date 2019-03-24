datafolder = ['C:\Users\User\Dropbox\Graduate Research Day 2019'...
    '\Optimization Data\Lattice Alt lat 30'];
walkerFolder = ['C:\Users\User\Dropbox\Graduate Research Day 2019'...
    '\Optimization Data\Walker lat 30'];
%% Choose Constellations
latEm = 30;
lonEm = 35;
hA = 900;
nSats = 63;
nPlanes = 9;

nSatsW = 63;
nPlanesW = 9;
Arch.nDays = 1;
Arch.nRepeats = 14; 
load ([datafolder '\OptParams.mat']);
load([datafolder '\LatticeExSol_Lat_' num2str(latEm) ...
    '_nSats_' num2str(nSats) '_hA_' num2str(hA) '.mat']);
%% Create Constellation
Arch.nSats = nSats;
Arch.nPlanes = nPlanes;
iPlane = ExSol.archMat(1,:) == nPlanes;
Arch.nAops = ExSol.archMat(2,iPlane);
Arch.nSatsPerAop = ExSol.archMat(3,iPlane);

Phase.nC1 = ExSol.phaseMat(1,iPlane);
Phase.nC2 = ExSol.phaseMat(2,iPlane);
Phase.nC3 = ExSol.phaseMat(3,iPlane);
ExSol.InitCon.raan1 = ExSol.InitCon.raan1 + lonEm;
LC = LatticeConstellation(Arch,Phase,ExSol.Orbit,ExSol.InitCon);
load([walkerFolder '\WalkerRgtExSol_Lat_' num2str(latEm)...
                '_T_' num2str(nSatsW) '.mat']);
ExSol.raan0 = ExSol.raan0 + lonEm;
[~,iOpt] = min(ExSol.intPdop(:,nPlanesW));
phasingFW = iOpt - 1;
WC = WalkerConstellation(nSatsW,nPlanesW,phasingFW,ExSol.inc,ExSol.alt,ExSol.raan0);
%% Propagate
Prop = Propagator(LC,PropParams.relTol,PropParams.absTol);
PropW = Propagator(WC,PropParams.relTol,PropParams.absTol);
[propTime, propState] = Prop.PropEciJ2(PropParams.timeVec);
[propTimeW, propStateW] = PropW.PropEciJ2(PropParams.timeVec);
% Evaluate PDOP
[pdop, satsIs] = TdoaPdopVec(propState,propTime,latEm,lonEm,0,PropParams.elevMin);
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
dLon = 10;
dLat = 10;
lats = max([5,latEm-dLat]):1:min([latEm+dLat,85]);
lons = [-dLon:1:dLon] + lonEm;
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
        pdop(pdop>100) = 1000;
        pdop(isnan(pdop)) = 1000;
        PDOP2(iLat) = trapz(propTime,pdop)/propTime(end);
        
        pdopW = TdoaPdopVec(propStateW,propTimeW,LAT(iLat,iLon),LON(iLat,iLon)...
            ,0,PropParams.elevMin);
        pdopW(pdopW>1000) = 1000;
        pdopW(isnan(pdopW)) = 1000;
        PDOPW2(iLat) = trapz(propTimeW,pdopW)/propTimeW(end);
    end
    PDOP(:,iLon) = PDOP2;
    PDOPW(:,iLon) = PDOPW2;
end
toc
%% Actually Plot
figure(2)
% subplot(2,1,1)

contourf(LON,LAT,PDOP,0:0.1:20,'LineColor','none')
% title(['Lattice Flower Constellation ' num2str(nSats) '/' num2str(nPlanes)])
hold on

colormap jet
shading interp
% c = colorbar;
% c.Label.String = 'Mean PDOP';
% c.FontSize = 10;
% c.Position = [0.6 0.1099 0.05 0.8154];
xlabel('Longitude')
ylabel('Latitude')
load coastlines
geoshow(coastlat,coastlon,'color','m','linewidth',1.5)
plot(lonEm,latEm,'yo','LineWidth',1.5)
axis equal
xlim([lonEm-dLon,lonEm+dLon])
xticks([lonEm-dLon:5:lonEm+dLon])
ylim([latEm-dLat,latEm+dLat])
yticks([latEm-dLat:5:latEm+dLat])
% grid on
hold off

% subplot(2,1,2)
figure(3)

contourf(LON,LAT,PDOPW,0:0.1:20,'LineColor','none')
% title(['Walker Constellation ' num2str(nSatsW) '/' num2str(nPlanesW)])
hold on

colormap jet
shading interp
% c = colorbar;
% c.Label.String = 'Mean PDOP';
% xlabel('Longitude')
% ylabel('Latitude')
load coastlines
geoshow(coastlat,coastlon,'color','m','linewidth',1.5)
plot(lonEm,latEm,'yo','LineWidth',1.5)
axis equal
xlim([lonEm-dLon,lonEm+dLon])
xticks([lonEm-dLon:5:lonEm+dLon])
ylim([latEm-dLat,latEm+dLat])
yticks([latEm-dLat:5:latEm+dLat])
% grid on
hold off