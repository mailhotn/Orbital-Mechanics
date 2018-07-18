
T = 54;
latGs = 20;
lonGs = 0;
time = 0:10:86164;
datafolder = 'C:\Users\User\Dropbox\Walker Optimization Data\Previous Optimization Runs\Walker GA Optimization RGT n20to80 j14to16';
%% Plot Nominal PDOP over Time
load([datafolder '\WalkerMeanRgtSolLat_' num2str(latGs) ...
                  '_T_' num2str(T) '.mat'])
alt = CalcRgtElement([],0,GaRgtSol.inc,GaRgtSol.jRepeats,1)-6378.137;
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
lat = 5:5:89;
lon = -90:10:90;
% lat = max([latGs-10,1]):min([latGs+10,89]);
% lon = -10:10;
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
contourf(LON,LAT,maxPdop,2000,'LineColor','none')
hold on
plot(lonGs,latGs,'y*')
colormap jet
shading interp
c = colorbar;
c.Label.String = 'Max PDOP';
xlabel('Longitude')
ylabel('Latitude')
title(['Max PDOP Map T = ' num2str(T) ' Ground Station Latitude =' ...
    num2str(latGs)])
hold off