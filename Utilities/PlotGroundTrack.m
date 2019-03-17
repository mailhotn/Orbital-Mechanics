function PlotGroundTrack( xEci, timeVec, gmst0 )
%plotGroundTrack Plots the ground tracks of a constellation
%   xEci is a Mx(6*N) matrix of ECI states
%   timeVec is a Mx1 or 1xM vector of time values
if nargin < 3
    gmst0 = 0;
end
primary = earth();
gmst = gmst0 + primary.we*timeVec;
nSats = size(xEci,2)/6;
tracks = nan(2*nSats,length(gmst));
for iTime = 1:length(gmst)
    xEcef = eci2ecef(reshape(xEci(iTime,:).',6,nSats),gmst(iTime));
    rEcef = xEcef(1:3,:)./sqrt(dot(xEcef(1:3,:),xEcef(1:3,:),1));
    lla   = ecef2lla(rEcef.'*primary.Re,0,primary.Re).';
    tracks(:,iTime) = reshape(lla(1:2,:),2*nSats,1);
end
gcf
hold on
load coastlines
geoshow(coastlat,coastlon)
for iSat = 1:nSats
    plot(wrapTo180(tracks(2*iSat,:)),tracks(2*iSat-1,:),'.');
end
axis equal
xlim([-180,180])
ylim([-90,90])
yticks(-90:30:90)
yticklabels({'90°S','60°S','30°S','0°','30°N','60°N','90°N'})
xticks(-180:60:180)
xticklabels({'180°W','120°W','60°W','0°','60°E','120°E','180°E'})
% xlabel('Longitude')
% ylabel('Latitude')
grid on
hold off
end

