function PlotPdopMap(Con,latEm,elevMin)
if nargin < 3
    elevMin = 5;
end
lonEm = 0;
timeVec = 0:100:86400;
Prop = Propagator(Con);
[propTime, propState] = Prop.PropEciJ2(timeVec);


%% Plot PDOP Map
dLon = 20;
dLat = 15;
lats = max([5,latEm-dLat]):1:min([latEm+dLat,85]);
lons = [-dLon:5:dLon] + lonEm;
[LON,LAT] = meshgrid(lons,lats);
PDOP = nan(size(LON));
for iLon = 1:length(lons)
    PDOP2 = zeros(length(lats),1);
    for iLat = 1:length(lats)
        pdop = TdoaPdopVec(propState,propTime,LAT(iLat,iLon),LON(iLat,iLon)...
            ,0,elevMin);
        pdop(pdop>1000) = 1000;
        pdop(isnan(pdop)) = 1000;
        PDOP2(iLat) = trapz(propTime,pdop)/propTime(end);
    end
    PDOP(:,iLon) = PDOP2;
end
%% Actually Plot
gcf

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
