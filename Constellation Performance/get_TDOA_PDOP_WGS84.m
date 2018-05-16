function [ PDOP ] = get_TDOA_PDOP_WGS84( X_GS, X_ECEF)
%get_TDOA_PDOP_sphere calculates PDOP of TDOA sats relative to Ground station
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
if size(X_ECEF,2) < 3
    PDOP = 20;
    return
end
e2 = 0.0818191908426215^2; % eccentricity squared
R  = X_ECEF(1:3,:);
dR = X_GS-R;
dR_hat = dR./sqrt(dot(dR,dR,1));

H = (dR_hat(:,2:end) - dR_hat(:,1)).';
PHI = [eye(2);
      -(1-e2)*X_GS(1)/X_GS(3),-(1-e2)*X_GS(2)/X_GS(3)];
G = PHI.'*(H.'*H)*PHI;
if cond(G) <= 1e10
    PDOP = sqrt(trace(PHI*inv(PHI.'*(H.'*H)*PHI)*PHI.')); %#ok<MINV>
else
    PDOP = 20;
end
end