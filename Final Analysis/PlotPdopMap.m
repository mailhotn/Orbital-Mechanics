function PlotPdopMap(Con,latEm,lonEm,elevMin,parallel,colorBar,prct)
if nargin < 3
    lonEm = 0;
    elevMin = 5;
    parallel = 0;
elseif nargin < 4
    elevMin = 5;
    parallel = 0;
elseif nargin < 5
    parallel = 0;
elseif nargin < 6
    colorBar = false;
end

timeVec = 0:100:86400;
Prop = Propagator(Con);
[propTime, propState] = Prop.PropEciJ2(timeVec);


%% Plot PDOP Map
dLon = 30;
dLat = 20;
lats = max([5,latEm-dLat]):1:min([latEm+dLat,85])+1;
lons = [-dLon:1:dLon+1] + lonEm;
[LON,LAT] = meshgrid(lons,lats);
PDOP = nan(size(LON));

if parallel
    parfor iLon = 1:length(lons)
        PDOP2 = zeros(length(lats),1);
        for iLat = 1:length(lats)
            pdop = TdoaPdopVec(propState,propTime,LAT(iLat,iLon),LON(iLat,iLon)...
                ,0,elevMin);
            pdop(pdop>1000) = 1000;
            pdop(isnan(pdop)) = 1000;
%             p95 = prctile(pdop,95);
%             pdop(pdop > p95) = nan;
%             PDOP2(iLat) = trapz(propTime(~isnan(pdop)),pdop(~isnan(pdop)))/propTime(end);
            PDOP2(iLat) = prctile(pdop,prct);
        end
        PDOP(:,iLon) = PDOP2;
    end
else
    for iLon = 1:length(lons)
        PDOP2 = zeros(length(lats),1);
        for iLat = 1:length(lats)
            pdop = TdoaPdopVec(propState,propTime,LAT(iLat,iLon),LON(iLat,iLon)...
                ,0,elevMin);
            pdop(pdop>1000) = 1000;
            pdop(isnan(pdop)) = 1000;
%             PDOP2(iLat) = trapz(propTime,pdop)/propTime(end);
            PDOP2(iLat) = prctile(pdop,prct);
        end
        PDOP(:,iLon) = PDOP2;
    end
end
PDOP(end,end) = 0;
PDOP(1,end) = 20;
%% Actually Plot
gcf

contourf(LON,LAT,PDOP,0:0.1:20,'LineColor','none')
% title(['Lattice Flower Constellation ' num2str(nSats) '/' num2str(nPlanes)])
hold on

colormap jet
shading interp
if colorBar
    c = colorbar;
    c.Label.String = "$p_{" + prct + "}\left(PDOP\right)$";
    c.Label.Interpreter = 'latex';
    c.Label.FontSize = 14;
end
% c = colorbar;
% c.Label.String = 'Mean PDOP';
% c.FontSize = 10;
% c.Position = [0.6 0.1099 0.05 0.8154];
% xlabel('Longitude')
% ylabel('Latitude')
load coastlines
geoshow(coastlat,coastlon,'color','m','linewidth',1.5);
plot(lonEm,latEm,'yo','LineWidth',1.5);
axis equal
xlim([lonEm-dLon,lonEm+dLon])
xticks([lonEm-dLon:10:lonEm+dLon])
ylim([latEm-dLat,latEm+dLat])
yticks([latEm-dLat:10:latEm+dLat])
% grid on
hold off
