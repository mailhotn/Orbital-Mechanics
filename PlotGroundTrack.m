function PlotGroundTrack( xEci, timeVec, gmst0 )
%plotGroundTrack Plots the ground tracks of a constellation
%   Detailed explanation goes here
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
    plot(wrapTo180(tracks(2*iSat,:)),tracks(2*iSat-1,:),'.','LineWidth',1.5)
end
axis equal
xlim([-180,180])
ylim([-90,90])
grid on
hold off
end

