function [rGsEst, pdop, exitFlag] = GeolocateSphere(rGs, xEcef, timeVar)
%GeolocateSphere estimates the position of the ground station
%   based on the satellite locations, the true ground station position and
%   the noise variance.
% -Based on a Spherical Earth model
% -Algorithm from:  K. C. Ho and Y. T. Chan, “Geolocation of a known
%  altitude object from TDOA and FDOA measurements,” 1997
% -Solution with fmincon, not Newton
%% Input Arguments
% * xGs     - 3x1 position vector of ground station in ECEF frame
% * xEcef   - 6xN matrix of satellite states in ECEF frame
% * timeVar - Variance of time measurement [sec^2]
% * elevMin - Minimum elevation angle (used for ambiguity resolution) [deg]
%% Exit Flag Values
% 1 - No ambiguity in solution
% 2 - Ambiguity resolved by measurement equation
% 3 - Ambiguity resolved by line of sight check
% 4 - No solution
% 5 - Ambiguity unresolved

%% Check #Sats
if size(xEcef,2) < 3
    rGsEst = nan(3,1);
    pdop = nan;
    exitFlag = 4;
    return
end
%% Constants
c = 3e5; % speed of light [km/s]
e2 = 0;
% e2 = 0.0818191908426215^2; % eccentricity squared (WGS84)
% conMat = diag([1,1,1/(1-e2)]);
%% True Ranges & Differences
rSats = xEcef(1:3,:);
nTdoas = size(xEcef,2) - 1;
range = sqrt(dot(rSats - rGs,rSats - rGs,1)).';
rdoaTrue = range(2:end)-range(1);
%% Add Gaussian Noise
Q = timeVar*(eye(nTdoas) + ones(nTdoas));
noise = mvnrnd(zeros(1,nTdoas),Q).';
rdoa = rdoaTrue + c*noise;

%% Estimation
options = optimoptions('fmincon','display','none');
primary = earth();
lb = [-inf,-inf,-inf,0];
ub = [primary.Re,primary.Re,primary.Re,inf];
x0 = [rSats(:,1);norm(rSats(:,1))];
x = fmincon(@(x)GeolocationCost(x,rSats,rdoa),x0,[],[],[],[],lb,ub,...
    @(x)GeolocationConstraints(x,rSats),options);

rGsEst = x(1:3);
pdop = 3;
exitFlag = 2;
end
function J = GeolocationCost(x,rSats,rdoa)
u = [x(1);x(2);x(3)];
r1 = x(4);
nTdoas = length(rdoa);
h = rdoa.^2 - dot(rSats(:,2:end),rSats(:,2:end),1).' +...
    dot(rSats(:,1),rSats(:,1),1).';
G1 = -2*(rSats(:,2:end) - rSats(:,1)).';
g2 = -2*rdoa;
W = eye(nTdoas) + ones(nTdoas);
J = (h - G1*u - g2*r1).'*W*(h - G1*u - g2*r1);
end
function [c,ceq] = GeolocationConstraints(x,rSats)
primary = earth();
u = [x(1);x(2);x(3)];
r1 = x(4);
ceq(1) = 2*rSats(:,1).'*u - rSats(:,1).'*rSats(:,1) - primary.Re^2 + r1^2;
ceq(2) = u.'*u - primary.Re^2;
c = [];
end