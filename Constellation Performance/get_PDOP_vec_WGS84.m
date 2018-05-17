function [ PDOP ] = get_PDOP_vec_WGS84(X, T, lat_gs, lon_gs, GMST0, e_min)
%get_PDOP_vec_sphere Calculates PDOP at each time step
% Based on WGS-84 Earth Model
%% Input Arguments
% * X      - Mx(6N) matrix of ECI states for N satellites at M time steps
% * T      - Mx1 vector of time steps
% * lat_gs - Ground Station latitude (deg)
% * lon_gs - Ground Station longitude (deg)
% * GMST0  - Initial GMST (deg)
% * e_min  - Minimum elevation for line of sight check
%
%%   ~~~~~~~~~~~~~ Algorithm ~~~~~~~~~~~~~
% 1. for loop over all time-steps, at each step:
%    a. Transform to ECEF
%    b. Find all satellites in line of sight 
%    c. Calculate PDOP at time step
% 2. Return Mx1 Vector of PDOP at each time step
% 
%% 
M = length(T);
N = size(X,2)/6;
w_e = 7.2921150e-5; % rad/sec
Re  = 6378.137; % km
GMST = wrapTo360(GMST0 + 180/pi*T*w_e);
GS = lla2ecef([lat_gs,lon_gs,0]).'/1000;
PDOP = 20*ones(M,1);
for ii = 1:M
    X_ECI = reshape(X(ii,:).',6,N);
    X_ECEF = eci2ecef(X_ECI,GMST(ii));
    X_IS = sats_in_sight(X_ECEF,GS,e_min,norm(GS));
    PDOP(ii) = get_TDOA_PDOP_WGS84(GS,X_IS);
end
end

