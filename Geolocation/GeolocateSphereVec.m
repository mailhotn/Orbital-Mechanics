function [rGsEst, estError, pdop, flagVec] = ...
    GeolocateSphereVec(xEci, time, latGs,lonGs, elevMin, timeVar)
%Geolocate3SatsVec Estimates the ground station position over time
%   Spherical Earth model
%% Input Arguments
% * xEci    - Mx(6N) matrix of ECI states for N satellites at M time steps
% * time    - Mx1 vector of time steps
% * latGs   - Ground Station latitude (deg)
% * lonGs   - Ground Station longitude (deg)
% * elevMin - Minimum elevation for line of sight check
% * timeVar - Variance of time measurement [sec^2]
primary = earth();
nTime = length(time);
nSats = size(xEci,2)/6;
gmst = wrapTo360(time*primary.we);
rGs = primary.Re*[cosd(lonGs)*cosd(latGs); sind(lonGs)*cosd(latGs); sind(latGs)];
flagVec = nan(nTime,1);
rGsEst = nan(6,nTime);
pdop = nan(1,nTime);
estError = inf(2,nTime);

for iTime = 1:nTime
    xEciNow = reshape(xEci(iTime,:).',6,nSats);
    xEcef = eci2ecef(xEciNow,gmst(iTime));
    xInSight = SatsInSight(xEcef,rGs,elevMin,primary.Re);
    [rGsNow, pdopNow, flagNow] = GeolocateSphere(rGs,xInSight,timeVar);
    flagVec(iTime) = flagNow;
    pdop(iTime) = pdopNow;
    if flagNow ~= 5
        rGsEst(1:3,iTime) = rGsNow;
        estError(1,iTime) = norm(rGsNow-rGs);
    else
        err1  = norm(rGsNow(1:3)-rGs);
        err2  = norm(rGsNow(4:6)-rGs);
        err12 = [err1,err2];
        [~,iReal] = min(err12);
        iFake = 3-iReal;
        rGsEst(:,iTime) = [rGsNow((iReal*3-2):(iReal*3));
                           rGsNow((iFake*3-2):(iFake*3))];
        estError(1,iTime) = err12(iReal);
        estError(2,iTime) = err12(iFake);
    end
end
end