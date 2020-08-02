function [minDist,minTime,minSats] = CalcMinDist(Con)
% CalcMinDist calculates the minimum distance between any 2 satellites in a
% constellation.
if Con.ecc == 0
    minDist = inf;
    oe0 = Con.InitialOeMean;
    for iSat = 1:(Con.nSats-1)
        for iSat2 = (iSat+1):Con.nSats
            delM = wrapTo360(oe0(6,iSat)+oe0(5,iSat)-oe0(6,iSat2)-oe0(5,iSat2));
            delO = wrapTo360(oe0(4,iSat)-oe0(4,iSat2));
            delF = delM -2*atand(-tand(delO/2)*cosd(Con.inc));
            cBeta = cosd(Con.inc)^2 + sind(Con.inc)^2*cosd(delO);
            rho = Con.sma*2*abs(sqrt((1+cBeta)/2)*sind(delF/2));
            if rho < minDist
                minDist = min([minDist, rho]);
                minSats = [iSat,iSat2];
                minTime = 0;
            end
        end
    end
else
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
end
