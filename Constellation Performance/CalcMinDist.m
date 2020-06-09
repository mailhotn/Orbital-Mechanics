function [minDist,minTime,minSats] = CalcMinDist(Con)
% CalcMinDist calculates the minimum distance between any 2 satellites in a
% constellation.
prop = Propagator(Con);

timeVec = 0:10:86164;
[~, xEci] = prop.PropEciJ2(timeVec);
xEci = reshape(xEci.',6,Con.nSats*length(timeVec));
rEci = xEci(1:3,:);
rEci = reshape(rEci,3*Con.nSats,length(timeVec));
minDist = inf;
for iSat = 1:Con.nSats
    satsList = iSat+1:Con.nSats;
    for iSat2 = satsList
        distVecs = rEci(3*iSat-2:3*iSat,:)-rEci(3*iSat2-2:3*iSat2,:);
        dist = sqrt(dot(distVecs,distVecs,1));
        if min(dist) < minDist
            [minDist,minTime] = min(dist);
            minSats = [iSat,iSat2];
        end
    end
end
