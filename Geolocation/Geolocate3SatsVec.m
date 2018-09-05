function [rGsEst, estError, pdop] = Geolocate3SatsVec(xEci, time, latGs,...
    lonGs, elevMin, timeVar)
%Geolocate3SatsVec Estimates the ground station position over time
%   Spherical Earth model, only 3 satellites taken into account
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
    [rGsNow, pdopNow, flagNow] = Geolocate3Sats(rGs,xInSight,timeVar,elevMin);
    flagVec(iTime) = flagNow;
    pdop(iTime) = pdopNow;
    if flagNow ~= 5
        rGsEst(1:3,iTime) = rGsNow;
        estError(1,iTime) = norm(rGsNow-rGs);
    else
        rGsEst(:,iTime) = rGsNow;
        estError(1,iTime) = norm(rGsNow(1:3)-rGs);
        estError(2,iTime) = norm(rGsNow(4:6)-rGs);
    end
end
end