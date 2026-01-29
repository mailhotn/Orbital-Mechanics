function dX = DynEci3B(P,t,X,ForceModel) 
% 3b acceleration from Vallado (8-34) fifth edition
% Normalized Parameters
mu = P.Con.primary.mu;
Re = P.Con.primary.Re;

% move stuff around
order = 6*P.Con.nSats;
X2 = reshape(X,[6,P.Con.nSats]);
Rv = [eye(3),zeros(3);
    zeros(3,6)]*X2;
% get vector of R magnitudes
%             r = sqrt(dot(Rv,Rv,1));
r = vecnorm(Rv);
r = reshape(repmat(r,6,1),order,1); % norm(R) in all 6 spaces
% move more stuff around
R2 = reshape(Rv,order,1);  % 6nx1 for dx
V2 = reshape([zeros(3,6);zeros(3),eye(3)]*X2,order,1);

% Third Bodies
n3 = length(P.Con.third);
f3B = zeros(order,1);
for i3 = 1:n3
    mu3 = P.Con.third{i3}.mu;
    % r3bVec = planetEphemeris(P.Con.epoch,P.Con.primary.name,P.Con.third{i3}.name).'; % Call PE each time
    % r3bVec = P.Con.x3AtEpoch(1:3,i3); % call PE at startup
    r3bVec = P.Con.third{i3}.PosJ2000(t); % Use approximate position
    r3b = norm(r3bVec); % scalar - 3b dist from primary
    R3b = repmat([r3bVec;zeros(3,1)],P.Con.nSats,1);  % 3b pos from primary; 6nx1 for dx
    R3bM = reshape(R3b,6,P.Con.nSats); % 6xn for dot
    RSat3 = R3b - R2; % 3b pos from sat 6nx1 for dx
    RSat3M = R3bM - Rv; % 6xn for dot
    rSat3 = vecnorm(([r3bVec;zeros(3,1)]-Rv)); % 3b dist from sat nx1
    rSat3 = reshape(repmat(rSat3,6,1),order,1); % 3b dist from sat repeated 6nx1 for dx
    rSatDotR3 = reshape(repmat(dot(Rv,R3bM,1),6,1),order,1); % sat pos dot 3b pos repeated 6nx1
    rSatDotRSat3 = reshape(repmat(dot(Rv,RSat3M,1),6,1),order,1); % sat pos 1 dot sat pos 3b repeated 6nx1

    Q = (r.^2 + 2*(rSatDotRSat3)).*(r3b^2 + r3b*rSat3 + rSat3.^2)./(r3b^3*rSat3.^3.*(r3b + rSat3));
    switch ForceModel
        case 'Taylor'
            f3B = -mu3./r3b.^3.*(R2 - 3*R3b.*rSatDotR3/r3b^2 - 15/2*(rSatDotR3/r3b^2).^2.*R3b); % approx (8-35)
        case 'Stable'
            f3B = f3B +  mu3*(Q.*RSat3 - R2/r3b^3); % (8-34)
        case 'Direct'
            f3B = f3B + mu3*(RSat3./rSat3.^3 -R3b/r3b^3); % (8-33)
    end
end
% equation of motion

dX = circshift(V2,-3) + circshift(-mu*R2./r.^3 + f3B,3);
