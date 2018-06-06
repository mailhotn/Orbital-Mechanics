T = 40;
latGs = 40;
time = 0:10:86164;
datafolder = 'C:\Users\User\Dropbox\Walker Optimization Data';
%% Plot Nominal PDOP over Time
load([datafolder '\WalkerMeanRgtSolLat_' num2str(latGs) ...
                  '_T_' num2str(T) '.mat'])
alt = CalcRgtElement([],0,GaRgtSol.inc,12,1)-6378.137;
WC = WalkerConstellation(GaRgtSol.nSatsT,GaRgtSol.nPlanesP,GaRgtSol.phasingF,...
    GaRgtSol.inc,alt);
Prop = Propagator(WC,1e-8,1e-9);
[propTime, propState] = Prop.PropEciJ2(time);
[pdop, satsIs] = TdoaPdopVec(propState,propTime,GaRgtSol.latGs,GaRgtSol.lonGs,...
                             GaRgtSol.gmst0, GaRgtSol.elevMin);
figure(1)
yyaxis left
plot(propTime/3600,pdop)
ylabel('PDOP')
yyaxis right
plot(propTime/3600,satsIs)
ylabel('Sats In Sight')
xlabel('Time [hr]')
grid
% plotyy(propTime(1:1000),pdop,propTime(1:1000),satsIs)

%% Plot Mean PDOP Map
lat = 2:2:40;
lon = -90:2:90;
[LON,LAT] = meshgrid(lon,lat);
maxPdop  = nan(size(LON));
tic
for ii = 1:size(LON,2)
    for jj = 1:size(LAT,1)
        pdopLocal = TdoaPdopVec(propState,propTime,LAT(jj,ii),LON(jj,ii),...
                             GaRgtSol.gmst0, GaRgtSol.elevMin);
        maxPdop(jj,ii) = max(pdopLocal);
    end
end
toc
figure(2)
contourf(LON,LAT,maxPdop,200,'LineColor','none')
% hold on
colormap jet
shading interp
c = colorbar;
c.Label.String = 'Mean PDOP';
xlabel('Longitude')
ylabel('Latitude')
title(['Mean PDOP Coverage T = ' num2str(T) ' Ground Station Latitude =' ...
    num2str(latGs)])