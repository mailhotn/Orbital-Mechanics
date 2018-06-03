function [ pdop ] = CalcTdoaPdop( xGs, xEcef)
%CalcTdoaPdop calculates PDOP of TDOA sats relative to Ground station
%   Based on WGS-84 Earth model
%% Input Arguments
%
% * X_GS   - 3x1 position vector of ground station in ECEF frame
% * X_ECEF - 6xN matrix of satellite states in ECEF frame
% * Re     - Primary radius (km)
%
%
% N - Number of satellites in view
%% Calculate PDOP
if size(xEcef,2) < 3
    pdop = 20;
    return
end
e2 = 0.0818191908426215^2; % eccentricity squared
R  = xEcef(1:3,:);
dR = xGs-R;
dR_hat = dR./sqrt(dot(dR,dR,1));

H = (dR_hat(:,2:end) - dR_hat(:,1)).';
PHI = [eye(2);
    -(1-e2)*xGs(1)/xGs(3),-(1-e2)*xGs(2)/xGs(3)];
G = PHI.'*(H.'*H)*PHI;
if cond(G) <= 1e10
    pdop = sqrt(trace(PHI*inv(PHI.'*(H.'*H)*PHI)*PHI.')); %#ok<MINV>
else
    pdop = 20;
end
if pdop > 20
    pdop = 20;
end
end