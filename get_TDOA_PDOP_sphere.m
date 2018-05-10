function [ PDOP ] = get_TDOA_PDOP_sphere( lat_gs, lon_gs, X_ECEF, Re)
%get_TDOA_PDOP_sphere calculates PDOP of TDOA sats relative to Ground station
%   Based on spherical Earth model
%% Input Arguments
% 
% * lat_gs - Ground Station latitude (deg)
% * lon_gs - Ground Station longitude (deg)
% * X_ECEF - 6xN matrix of satellite states in ECEF frame
%
%% Calculate PDOP
GS = Re*[cosd(lon_gs)*cosd(lat_gs);
         sind(lon_gs)*cosd(lat_gs);
         sind(lat_gs)];
R  = X_ECEF(1:3,:);
dR = GS-R;
dR_hat = dR./sqrt(dot(dR,dR,1));

H = [(dR_hat(:,2:end) - dR_hat(:,1)).'; R.'/Re];

PDOP = sqrt(trace(inv(H.'*H)));
end

