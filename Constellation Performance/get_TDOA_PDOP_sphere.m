function [ PDOP ] = get_TDOA_PDOP_sphere( X_GS, X_ECEF)
%get_TDOA_PDOP_sphere calculates PDOP of TDOA sats relative to Ground station
%   Based on spherical Earth model
%% Input Arguments
% 
% * X_GS   - 3x1 position vector of ground station in ECEF frame
% * X_ECEF - 6xN matrix of satellite states in ECEF frame
% * Re     - Primary radius (km)
% N - Number of satellites in view
%% Calculate PDOP
if size(X_ECEF,2) < 3
    PDOP = 20;
    return
end
R  = X_ECEF(1:3,:);
dR = X_GS-R;
dR_hat = dR./sqrt(dot(dR,dR,1));
H = [(dR_hat(:,2:end) - dR_hat(:,1)).'; X_GS.'];
G = H.'*H;
if cond(G) <= 1e10
    PDOP = sqrt(trace(inv(G)));
else
    PDOP = 20;
end
end