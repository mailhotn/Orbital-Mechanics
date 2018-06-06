function [ pdop , nSatsInSight] = TdoaPdopVec(X, T, latGs, lonGs, gmst0, elevMin)
%TdoaPdopVec Calculates PDOP at each time step
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
Earth = earth();
M = length(T);
N = size(X,2)/6;
w_e = Earth.we;
GMST = wrapTo360(gmst0 + T*w_e);
GS = lla2ecef([latGs,lonGs,0]).'/1000;
pdop = nan*ones(M,1);
nSatsInSight = zeros(M,1);
for iTime = 1:M
    X_ECI = reshape(X(iTime,:).',6,N);
    X_ECEF = eci2ecef(X_ECI,GMST(iTime));
    X_IS = SatsInSight(X_ECEF,GS,elevMin,norm(GS));
    nSatsInSight(iTime) = size(X_IS,2);
    pdop(iTime) = CalcTdoaPdop(GS,X_IS);
end
end

