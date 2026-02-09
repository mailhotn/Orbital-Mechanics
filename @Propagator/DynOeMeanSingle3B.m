function dX = DynOeMeanSingle3B(P,t,X)
% Singly averaged dynamics with third bodies
% Based on T. Nie and P. Gurfil, “Long-term evolution of 
% orbital inclination due to third-body inclination,”  Jan. 2021
oe  = reshape(X,6,P.Con.nSats);
mu = P.Con.primary.mu;
% Element vectors (angles already in radians)
sma = oe(1,:);
ecc = oe(2,:);
inc = oe(3,:);
ran = oe(4,:);
aop = oe(5,:);
man = oe(6,:);
nMo = sqrt(mu./sma.^3);
eta = sqrt(1-ecc.^2);

%% Calculate element rate for each third body

n3 = length(P.Con.third);
dSma = zeros(1,P.Con.nSats);
dEcc = zeros(1,P.Con.nSats);
dInc = zeros(1,P.Con.nSats);
dRan = zeros(1,P.Con.nSats);
dAop = zeros(1,P.Con.nSats);
dMan = zeros(1,P.Con.nSats);

for j3 = 1:n3
    mu3 = P.Con.third{j3}.mu;
    % [r3b,v3b] = planetEphemeris(P.Con.epoch,P.Con.primary.name,P.Con.third{j3}.name);
    % oe3 = eci2oe3b(r3b.',v3b.',P.Con.primary,P.Con.third{j3});
    % x3b = P.Con.x3AtEpoch;
    [r3bMat, v3bMat] = P.Con.third{j3}.PosJ2000(t);
    % r3b = x3b(1:3,j3);
    % v3b = x3b(4:6,j3);
    oe3 = eci2oe3b(r3bMat,v3bMat,P.Con.primary,P.Con.third{j3});

    r3 = vecnorm(r3bMat,2,1);
    a3 = oe3(1,:);
    e3 = oe3(2,:);
    i3 = oe3(3,:);
    O3 = oe3(4,:);
    u3 = oe3(5,:)+oe3(6,:);

    th = ran-O3;

    alpha = (cos(aop).*cos(th) - sin(aop).*cos(inc).*sin(th)).*cos(u3) + ...
        (sin(aop).*sin(inc).*sin(i3) + cos(aop).*cos(i3).*sin(th) + sin(aop).*cos(inc).*cos(i3).*cos(th)).*sin(u3);

    beta = (-sin(aop).*cos(th) - cos(aop).*cos(inc).*sin(th)).*cos(u3) + ...
        (cos(aop).*sin(inc).*sin(i3) - sin(aop).*cos(i3).*sin(th) + cos(aop).*cos(inc).*cos(i3).*cos(th)).*sin(u3);

    dadi = sin(aop).*((-sin(u3).*cos(i3).*cos(th) + sin(th).*cos(u3)).*sin(inc) + cos(inc).*sin(i3).*sin(u3));
    dbdi = cos(aop).*((-sin(u3).*cos(i3).*cos(th) + sin(th).*cos(u3)).*sin(inc) + cos(inc).*sin(i3).*sin(u3));

    dado = ((cos(inc).*cos(i3).*cos(th) + sin(inc).*sin(i3)).*sin(u3) - sin(th).*cos(inc).*cos(u3)).*cos(aop) - sin(aop).*(sin(th).*sin(u3).*cos(i3) + cos(th).*cos(u3));
    dbdo = ((-cos(inc).*cos(i3).*cos(th) - sin(inc).*sin(i3)).*sin(u3) + sin(th).*cos(inc).*cos(u3)).*sin(aop) - cos(aop).*(sin(th).*sin(u3).*cos(i3) + cos(th).*cos(u3));

    dadth = (sin(u3).*cos(i3).*cos(th) - sin(th).*cos(u3)).*cos(aop) - cos(inc).*sin(aop).*(sin(th).*sin(u3).*cos(i3) + cos(th).*cos(u3));
    dbdth = -cos(inc).*(sin(th).*sin(u3).*cos(i3) + cos(th).*cos(u3)).*cos(aop) - sin(aop).*(sin(u3).*cos(i3).*cos(th) - sin(th).*cos(u3));

    k3 = mu3/4./r3.^3./nMo;
    Rm = k3.*(3*alpha.^2.*(4*ecc.^2+1)-3*beta.^2.*(ecc.^2-1)-3*ecc.^2-2);
    dRda = Rm./sma*2;
    dRde = k3.*(3*alpha.^2.*(8*ecc)-3*beta.^2.*(2*ecc)-6*ecc);
    dRdi = k3.*(6*alpha.*dadi.*(4*ecc.^2+1)-6*beta.*dbdi.*(ecc.^2-1));
    dRdO = k3.*(6*alpha.*dadth.*(4*ecc.^2+1)-6*beta.*dbdth.*(ecc.^2-1));
    dRdo = k3.*(6*alpha.*dado.*(4*ecc.^2+1)-6*beta.*dbdo.*(ecc.^2-1));
    dSma = dSma + zeros(1,P.Con.nSats);
    dEcc = dEcc - eta./ecc.*dRdo;
    dInc = dInc + 1./(eta.*sin(inc)).*(cos(inc).*dRdo - dRdO);
    dRan = dRan + 1./(eta.*sin(inc)).*dRdi;
    dAop = dAop + eta./ecc.*dRde -cot(inc)./eta.*dRdi;
    dMan = dMan - eta./ecc.*dRde - 2*sma.*dRda;
    
end

dOe = reshape([dSma;dEcc;dInc;dRan;dAop;dMan + nMo],6*P.Con.nSats,1);
% Equations of Motion
dX = dOe;