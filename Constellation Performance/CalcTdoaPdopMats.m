function [ measMat, conMat, Qinv ] = CalcTdoaPdopMats( xGs, xEcef)
%CalcTdoaPdop calculates PDOP of TDOA sats relative to Ground station
%   Based on WGS-84 Earth model
%% Input Arguments
%
% * X_GS   - 3x1 position vector of ground station in ECEF frame
% * X_ECEF - 6xN matrix of satellite states in ECEF frame
%
%
% N - Number of satellites in view
%% Calculate PDOP
Qinv = inv(eye(size(xEcef,2)-1) + ones(size(xEcef,2)-1));
e2 = 0.0818191908426215^2; % eccentricity squared
R  = xEcef(1:3,:);
dR = xGs-R;
dR_hat = dR./sqrt(dot(dR,dR,1));

measMat = (dR_hat(:,2:end) - dR_hat(:,1)).';
conMat = [eye(2);
    -(1-e2)*xGs(1)/xGs(3),-(1-e2)*xGs(2)/xGs(3)];


%% Other Version - Numerically unstable, but not singular for latGs = 0
% condTol = 1e15;
% %% Calculate PDOP
% if size(xEcef,2) < 3
%     pdop = nan;
%     return
% end
% Qinv = inv(eye(size(xEcef,2)-1) + ones(size(xEcef,2)-1));
% 
% R  = xEcef(1:3,:);
% dR = xGs-R;
% dR_hat = dR./sqrt(dot(dR,dR,1));
% H = (dR_hat(:,2:end) - dR_hat(:,1)).';
% 
% fish = (H.'*Qinv*H);
% 
% if cond(fish) <= condTol
%     e2 = 0.0818191908426215^2; % eccentricity squared
%     w = [xGs(1:2);xGs(3)/(1-e2)];
%     
%     Xinv = inv(fish);
%     G = Xinv - (Xinv*(w*w.')*Xinv)./(w.'*Xinv*w); %#ok<*MINV>
%     pdop = sqrt(trace(G));
% else
%     pdop = nan;
% end