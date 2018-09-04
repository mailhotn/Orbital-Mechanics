function [rGsEst, pdop, exitFlag] = Geolocate3Sats(rGs, xEcef, timeVar,...
                                                     elevMin)
%Geolocate3Sats estimates the position of the ground station
%   based on the satellite locations, the true ground station position and
%   the noise variance.
% -This algorithm only takes into account 3 satellites
% -Based on a Spherical Earth model
% -Algorithm from:  K. C. Ho and Y. T. Chan, “Geolocation of a known
%  altitude object from TDOA and FDOA measurements,” 1997
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
conMat = diag([1,1,1/(1-e2)]);
%% True Ranges & Differences
rSats = xEcef(1:3,1:3);
range = sqrt(dot(rSats - rGs,rSats - rGs,1)).';
rdoaTrue = range(2:3)-range(1);
%% Add Gaussian Noise
Q = timeVar*[2 1;
             1 2];
noise = mvnrnd([0;0],Q).';
rdoa = rdoaTrue + c*noise;
%% Calulate PDOP
H = -[(rSats(:,2)-rGs).'/range(2) - (rSats(:,1)-rGs).'/range(1);
      (rSats(:,3)-rGs).'/range(3) - (rSats(:,1)-rGs).'/range(1);
                                                          rGs.'];
pdop = sqrt(trace(H\[Q zeros(2,1);zeros(1,3)]/H.'/timeVar));
%% Estimation
G = -2*[              rSats(:,1).';
         rSats(:,2).'-rSats(:,1).';
         rSats(:,3).'-rSats(:,1).'];
h = [                      -rGs.'*rGs-rSats(:,1).'*rSats(:,1)         0 1;
    rdoa(1)^2-rSats(:,2).'*rSats(:,2)+rSats(:,1).'*rSats(:,1) 2*rdoa(1) 0;
    rdoa(2)^2-rSats(:,3).'*rSats(:,3)+rSats(:,1).'*rSats(:,1) 2*rdoa(2) 0];
M = (h.'/G.')*conMat*(G\h);
solPoly = [M(3,3), M(2,3)+M(3,2), M(1,3)+M(2,2)+M(3,1),...
           M(1,2)+M(2,1) , M(1,1)-rGs.'*rGs].';
r1Sol = roots(solPoly).';
% Find actual solution
isPos = real(r1Sol) > 0;
isReal = imag(r1Sol) == 0;
r1Sol = r1Sol(isReal & isPos);
polyVec = [ones(1,length(r1Sol)); r1Sol; r1Sol.^2];

xSol = G\(h*polyVec);
if numel(xSol(~isnan(xSol))) == 3
    exitFlag = 1;
    rGsEst = xSol(~isnan(xSol));
    return
end
for iSol = 1:size(xSol,2)
    rangeEst = sqrt(dot(rSats - xSol(:,iSol),rSats - xSol(:,iSol),1)).';
    rdoaEst = rangeEst(2:3) - rangeEst(1);
    if norm(rdoaEst - rdoa) > 1e-8
        xSol(:,iSol) = nan(3,1);
    end
end
if numel(xSol(~isnan(xSol))) == 3
    exitFlag = 2;
    rGsEst = xSol(~isnan(xSol));
    return
end
% Resolve Ambiguity with LOS
for iSol = 1:size(xSol,2)
    satsInLos = SatsInSight(rSats,xSol(:,iSol),elevMin,norm(xSol(:,iSol)));
    if size(satsInLos,2) < 3
        xSol(:,iSol) = nan(3,1);
    end
end
if numel(xSol(~isnan(xSol))) == 3
    exitFlag = 3;
    rGsEst = xSol(~isnan(xSol));
    return
elseif numel(xSol(~isnan(xSol))) == 0
    exitFlag = 4;
    rGsEst = nan(3,1);
    return
else
    exitFlag = 5;
    rGsEst = xSol(~isnan(xSol));
end