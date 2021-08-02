e = 0.99;
M = 0:0.001:2*pi;
tol = 1e-14;
%% Calculate N-R, Battin
Etrue = M2E(M,e,tol);

cEtrue = cos(2*Etrue);

% cEBatt = CosEBatt(M,e,tol);

cEAlt = Cos2E(M,e,tol);
error = norm(cEAlt-cEtrue)

%% Plot Stuff
figure(1) % M & E
plot(M,Etrue,'linewidth',2)
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
yticks(0:pi/2:2*pi)
yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
axis([0 2*pi 0 2*pi])

figure(2)
plot(M,cEtrue,M,cEAlt,'--')
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
axis([0 2*pi -1 1])


function cE = CosEBatt(M,e,tol)
if nargin < 3
    tol = 1e-14;
end
cE = nan(1,length(M));
for iM = 1:length(M)
    cE(iM) = -e/2;
    dE = inf;
    k = 1;
    while abs(dE) > tol
        dE = (1/k)*(besselj(k-1,k*e) - besselj(k+1,k*e))*cos(k*M(iM));
        cE(iM) = cE(iM) + dE;
        k = k+1;
    end
end
end

function cE = Cos2E(M,e,tol)
if nargin < 3
    tol = 1e-14;
end
cE = nan(1,length(M));
for iM = 1:length(M)
    cE(iM) = 0;
    dE = inf;
    k = 1;
    while abs(dE) > tol
        dE = (2/k)*(besselj(k-2,k*e) - besselj(k+2,k*e))*cos(k*M(iM));
        cE(iM) = cE(iM) + dE;
        k = k+1;
    end
end
end