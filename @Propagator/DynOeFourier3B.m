function [freq0,lpeSpec] = DynOeFourier3B(P,t,icM,kMax)
%% Handle Input
J2 = P.Con.primary.J2;
Re = P.Con.primary.Re;
mu = P.Con.primary.mu;

% handle elements vector
nT = length(t);

% handle elements vector
sma = icM(1);
e = icM(2);
inc = icM(3);
ran = icM(4);
aop = icM(5);

nMo = sqrt(mu/sma^3);
p = sma*(1-e^2);
eta = sqrt(1-e^2);


lpeSpec = zeros(12,kMax+1);

%% Calculate Fourier coefficients independant of 3B
aVec = [e, -4*e^2-2, 4*e^3+11*e, -16*e^2-4, 4*e^3+11*e, -4*e^2-2, e];
bVec = [e, -2*e^2-2, 5*e, 0, -5*e, 2*e^2+2, -e];
cVec = [e^3, -6*e^2, 3*e^3+12*e, -12*e^2-8, 3*e^3+12*e, -6*e^2, e^3];

daVec = [1, -8*e, 12*e^2+11, -32*e, 12*e^2+11, -8*e, 1];
dbVec = [1, -4*e, 5, 0, -5, 4*e, -1];
dcVec = [3*e^2, -12*e, 9*e^2+12, -24*e, 9*e^2+12, -12*e, 3*e^2];

AkM = nan(1,kMax);
AkdeM = nan(1,kMax);

BkM = nan(1,kMax);
BkdeM = nan(1,kMax);

CkM = nan(1,kMax);
CkdeM = nan(1,kMax);

for k = 1:kMax
    n = (-k-3):(3-k);
    jVec = (besselj(n,-k*e)).';
    djVec = -0.5*k*(besselj(n-1,-k*e) - besselj(n+1,-k*e)).';

    AkM(k) = -0.25*aVec*jVec;
    BkM(k) = -eta/4*bVec*jVec;
    CkM(k) = -0.25*cVec*jVec;

    AkdeM(k) = -0.25*(daVec*jVec + aVec*djVec);
    BkdeM(k) = e/4/eta*bVec*jVec-eta/4*(dbVec*jVec + bVec*djVec);
    CkdeM(k) = -0.25*(dcVec*jVec + cVec*djVec);
end
% Freq 0 Ak, Bk, Ck
AkM = [(4*e^2+1)/2, AkM];
AkdeM = [4*e, AkdeM];
BkM = [0, BkM];
BkdeM = [0, BkdeM];
CkM = [(3*e^2+2)/2, CkM];
CkdeM = [3*e,CkdeM];

%% Calculate spectrum for each third body
n3 = length(P.Con.third);
for j3 = 1:n3
    mu3 = P.Con.third{j3}.mu;
    [r3b,v3b] = planetEphemeris(P.Con.epoch,P.Con.primary.name,P.Con.third{j3}.name);
    oe3 = eci2oe3b(r3b.',v3b.',P.Con.primary,P.Con.third{j3});
    
    r3 = norm(r3b);
    a3 = oe3(1);
    e3 = oe3(2);
    i3 = oe3(3);
    O3 = oe3(4);
    u3 = oe3(5)+oe3(6);

    th = ran-O3;

    alpha = (cos(aop)*cos(th) - sin(aop)*cos(inc)*sin(th))*cos(u3) + ...
        (sin(aop)*sin(inc)*sin(i3) + cos(aop)*cos(i3)*sin(th) + sin(aop)*cos(inc)*cos(i3)*cos(th))*sin(u3);

    beta = (-sin(aop)*cos(th) - cos(aop)*cos(inc)*sin(th))*cos(u3) + ...
        (cos(aop)*sin(inc)*sin(i3) - sin(aop)*cos(i3)*sin(th) + cos(aop)*cos(inc)*cos(i3)*cos(th))*sin(u3);

    dadi = sin(aop)*((-sin(u3)*cos(i3)*cos(th) + sin(th)*cos(u3))*sin(inc) + cos(inc)*sin(i3)*sin(u3));
    dbdi = cos(aop)*((-sin(u3)*cos(i3)*cos(th) + sin(th)*cos(u3))*sin(inc) + cos(inc)*sin(i3)*sin(u3));

    dado = ((cos(inc)*cos(i3)*cos(th) + sin(inc)*sin(i3))*sin(u3) - sin(th)*cos(inc)*cos(u3))*cos(aop) - sin(aop)*(sin(th)*sin(u3)*cos(i3) + cos(th)*cos(u3));
    dbdo = ((-cos(inc)*cos(i3)*cos(th) - sin(inc)*sin(i3))*sin(u3) + sin(th)*cos(inc)*cos(u3))*sin(aop) - cos(aop)*(sin(th)*sin(u3)*cos(i3) + cos(th)*cos(u3));

    dadth = (sin(u3)*cos(i3)*cos(th) - sin(th)*cos(u3))*cos(aop) - cos(inc)*sin(aop)*(sin(th)*sin(u3)*cos(i3) + cos(th)*cos(u3));
    dbdth = -cos(inc)*(sin(th)*sin(u3)*cos(i3) + cos(th)*cos(u3))*cos(aop) - sin(aop)*(sin(u3)*cos(i3)*cos(th) - sin(th)*cos(u3));



    % coefficients of series & derivatives

    c1 = 3*(alpha^2-beta^2);
    s1 = 6*alpha*beta;
    c0 = (3*beta^2-1);

    dc1di = 6*(alpha*dadi - beta*dbdi);
    ds1di = 6*(alpha*dbdi + beta*dadi);
    dc0di = 6*beta*dbdi;

    dc1do = 6*(alpha*dado - beta*dbdo);
    ds1do = 6*(alpha*dbdo + beta*dado);
    dc0do = 6*beta*dbdo;

    dc1dO = 6*(alpha*dadth - beta*dbdth);
    ds1dO = 6*(alpha*dbdth + beta*dadth);
    dc0dO = 6*beta*dbdth;

    k3 = P.Con.third{j3}.mass/(P.Con.primary.mass+P.Con.third{j3}.mass)*P.Con.third{j3}.nMo^2/2*a3^3/r3^3;

    % Common factors
    Rsma = 2*k3*sma/nMo;
    Recc = k3*eta/nMo/e;
    Rinc = k3/nMo/eta/sin(inc);
    Rran = k3/nMo/eta/sin(inc);
    Raop = k3/nMo;
    Rman = -k3/nMo;
    % Calculate Spectrum of Elements
    k = 0:kMax;

    lpeSpec(1:2,:) = lpeSpec(1:2,:) + ...
        Rsma*[s1*(BkM.*k);
        -(c1*AkM.*k+c0*CkM.*k)]; 
    lpeSpec(3:4,:) = lpeSpec(3:4,:) + ...
        Recc*[eta*s1*(BkM.*k) - (dc1do*AkM+dc0do*CkM);
        -eta*(c1*AkM.*k+c0*CkM.*k) - ds1do*BkM];
    lpeSpec(5:6,:) = lpeSpec(5:6,:) + ...
        Rinc*[cos(inc)*(dc1do*AkM + dc0do*CkM) - (dc1dO*AkM + dc0dO*CkM);
        cos(inc)*ds1do*BkM - ds1dO*BkM];
    lpeSpec(7:8,:) = lpeSpec(7:8,:) + ...
        Rran*[dc1di*AkM + dc0di*CkM;
        ds1di*BkM];
    lpeSpec(9:10,:) = lpeSpec(9:10,:) + ...
        Raop*[eta/e*(c1*AkdeM + c0*CkdeM) - cos(inc)/sin(inc)/eta*(dc1di*AkM + dc0di*CkM);
        eta/e*(s1*BkdeM) - cos(inc)/sin(inc)/eta*ds1di*BkM];
    lpeSpec(11:12,:) = lpeSpec(11:12,:) + ...
        Rman*[eta^2/e*(c1*AkdeM + c0*CkdeM) + 4*(c1*AkM + c0*CkM);
        eta^2/e*(s1*BkdeM) + 4*s1*BkM];
end
freq0 = lpeSpec(1:2:11,1);
lpeSpec = lpeSpec(:,2:end);

