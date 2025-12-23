function dX = DynEci3B(P,t,X,r3bVec) %#ok<INUSL>
% 3b acceleration from Vallado (8-36)
% Normalized Parameters
mu = P.Con.primary.mu;
Re = P.Con.primary.Re;
mu3 = P.Con.third.mu;

% move stuff around
order = 6*P.Con.nSats;
X2 = reshape(X,[6,P.Con.nSats]);
Rv = [eye(3),zeros(3);
    zeros(3,6)]*X2;
% get vector of R magnitudes
%             r = sqrt(dot(Rv,Rv,1));
r = vecnorm(Rv);
r = reshape(repmat(r,6,1),order,1);
% move more stuff around
R2 = reshape(Rv,order,1);
V2 = reshape([zeros(3,6);zeros(3),eye(3)]*X2,order,1);

% Third Body
r3b = norm(r3bVec);
R3b = repmat([r3bVec;zeros(3,1)],P.Con.nSats,1);
R3bM = reshape(R3b,6,P.Con.nSats);
RSat3 = R3b - R2;
rSat3 = vecnorm(([r3bVec;zeros(3,1)]-Rv));
rSat3 = reshape(repmat(rSat3,6,1),order,1);
rSatDotR3 = reshape(repmat(dot(Rv,R3bM,1),6,1),order,1);

f3B = -mu3./r3b.^3.*(R2 - 3*R3b.*rSatDotR3/r3b^2 - 15/2*(rSatDotR3/r3b^2).^2.*R3b);
% equation of motion

dX = circshift(V2,-3) + circshift(-mu*R2./r.^3 + f3B,3);
end