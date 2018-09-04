function [xIs] = SatsInSight( xEcef, xGs, elevMin, Re )
%SatsInSight Gets satellites in line of sight to Ground Station
%   Algorithm adapted from "Orbit and Constellation Design and Management"
%   by James R Wertz 2001, page 423.
%% Input Arguments 
% 
% * X_ECEF - 6xN Matrix of satellites states in ECEF frame
% * X_GS   - 3x1 position vector of ground station in ECEF frame
% * e_min  - minimum acceptable elevation angle (deg)
% * Re     - >>>>> LOCAL <<<<< Earth radius at ground station 

%%
N = size(xEcef,2);
S = xEcef(1:3,:);
P = xGs;
O = P - S;
O_hat = O./sqrt(dot(O,O,1));
D = min([-dot(S,O_hat,1) + sqrt(dot(S,O_hat,1).^2 - dot(S,S,1) + Re^2);
         -dot(S,O_hat,1) - sqrt(dot(S,O_hat,1).^2 - dot(S,S,1) + Re^2)],[],1);
diff = abs(sqrt(dot(O,O,1))-D);

P_hat = P./sqrt(dot(P,P,1));
P_hat = repmat(P_hat,1,N);

elev   = acosd(dot(P_hat,O_hat,1))-90;

sight = ((diff<1e-6).*(elev>=elevMin));

xIs = xEcef(:,sight==1);
end