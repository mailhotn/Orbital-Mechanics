function [freq0,lpeSpec] = DynOeFourierSimplified(P,t,icM,kMax)
%% Handle Input
J2 = P.Con.primary.J2;
Re = P.Con.primary.Re;
mu = P.Con.primary.mu;

% handle elements vector
% nT = length(t);

% Propagator now gives Mean elements
% icM = osc2meNum(icOsc); % change to numerical mean
% icM(3:end) = icM(3:end)*pi/180;
% icOsc(3:end) = icOsc(3:end)*pi/180;
% handle elements vector
a = icM(1);
e = icM(2);
i = icM(3);

nMo = sqrt(mu/a^3);
% p = a*(1-e^2);
aop = icM(5);

b = (1-sqrt(1-e^2))/e;


%% constant vectors
m = (0:4).';


a2 = [1,-4*e,4*e^2+2,-4*e,1].';

b1 = [-1,2*e,0,-2*e,1].';

da2de = [0,-4,8*e,-4,0].';
db1de = [0,2,0,-2,0].';

c1 = -3*sin(i)^2*cos(2*aop);
s1 = 3*sin(i)^2*sin(2*aop);
c0 = 3*sin(i)^2*cos(aop)^2-1;

dc1di_si = -6*cos(i)*cos(2*aop);
dc1do = 6*sin(i)^2*sin(2*aop);
dc1do_si = 6*sin(i)*sin(2*aop);

ds1di_si = 6*cos(i)*sin(2*aop);
ds1do = 6*sin(i)^2*cos(2*aop);
ds1do_si = 6*sin(i)*cos(2*aop);

dc0di_si = 6*cos(i)*cos(aop)^2;
dc0do = -6*sin(i)^2*cos(aop)*sin(aop);
dc0do_si = -6*sin(i)*cos(aop)*sin(aop);

AkM = nan(1,kMax);
Ak_eM = nan(1,kMax);
Akde_eM = nan(1,kMax);

BkM = nan(1,kMax);
Bk_eM = nan(1,kMax);
Bkde_eM = nan(1,kMax);

CkM = nan(1,kMax);
Ck_eM = nan(1,kMax);
Ckde_eM = nan(1,kMax);


k = 1;

while k <= kMax
    n = 0;

    N = m + n + k + 1;
    Nc = n + k + 1;

    g3 = abs(Nc).*b.^abs(Nc-1) + e/sqrt(1-e^2)*b.^abs(Nc);

    g5 = abs(N-2).*abs(N-1).*abs(N)/6.*b.^abs(N-3) + ...
        e*abs(N-1).*abs(N)/sqrt(1-e^2).*b.^abs(N-2) + ...
        5*e^2*abs(N)/2/(1-e^2).*b.^abs(N-1) + ...
        5*e^3/2/(1-e^2)^(3/2).*b.^abs(N);

    dg3deXe = abs(Nc).*(abs(Nc-1).*b.^(abs(Nc-1)) + ...
        e/sqrt(1-e^2)*b.^(abs(Nc)))/sqrt(1-e^2) + ...
        e*b.^abs(Nc)/(1-e^2)^(3/2);

    dg5deXe = abs(N).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(N) + ...
        abs(N-1).*(5*e^2/2/(1-e^2)*b.^abs(N-1) + ...
        abs(N-2).*(e/sqrt(1-e^2)*b.^abs(N-2) + ...
        abs(N-3)/6.*b.^abs(N-3))))/sqrt(1-e^2) + ...
        abs(N-1).*abs(N)*e/(1-e^2)^(3/2).*b.^abs(N-2) +...
        5*abs(N)*e^2/(1-e^2)^2.*b.^abs(N-1) +...
        15*e^3/2/(1-e^2)^(5/2)*b.^abs(N);

    Jn = besselj(n,-k*e);
    % Jn/e nonsingular
    if n~=0
        Jn_e = -k/2/n*(besselj(n+1,-k*e)+besselj(n-1,-k*e));
    else
        Jn_e = Jn/e;
    end
    % dJnde/e nonsingular
    if abs(n)~=1
        dJnde_e = k^2/4/(n^2-1)*(2*Jn + (n+1)*besselj(n-2,-k*e) - (n-1)*besselj(n+2,-k*e));
    else
        dJnde_e = -0.5*k/e*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
    end
    % Jn/e^2 nonsingular
    if abs(n)>1
        Jn_e2 = k^2/4/n/(n^2-1)*(2*n*Jn + (n+1)*besselj(n-2,-k*e) + (n-1)*besselj(n+2,-k*e));
    else
        Jn_e2 = Jn/e^2;
    end

    Ak = Jn*a2.'*g5;
    Ak_e = Jn_e*a2.'*g5;
    Akde_e = dJnde_e*a2.'*g5 + Jn_e*da2de.'*g5 + Jn_e2*a2.'*dg5deXe;

    Bk = Jn*b1.'*g5;
    Bk_e = Jn_e*b1.'*g5;
    Bkde_e = dJnde_e*b1.'*g5 + Jn_e*db1de.'*g5 + Jn_e2*b1.'*dg5deXe;

    Ck = Jn*g3;
    Ck_e = Jn_e*g3;
    Ckde_e = dJnde_e*g3  + Jn_e2*dg3deXe;


    n = 1;
    while n <= k + 3
        % positive n
        N = m + n + k + 1;
        Nc = n + k + 1;

        g3 = abs(Nc).*b.^abs(Nc-1) + e/sqrt(1-e^2)*b.^abs(Nc);

        g5 = abs(N-2).*abs(N-1).*abs(N)/6.*b.^abs(N-3) + ...
            e*abs(N-1).*abs(N)/sqrt(1-e^2).*b.^abs(N-2) + ...
            5*e^2*abs(N)/2/(1-e^2).*b.^abs(N-1) + ...
            5*e^3/2/(1-e^2)^(3/2).*b.^abs(N);

        dg3deXe = abs(Nc).*(abs(Nc-1).*b.^(abs(Nc-1)) + ...
            e/sqrt(1-e^2)*b.^(abs(Nc)))/sqrt(1-e^2) + ...
            e*b.^abs(Nc)/(1-e^2)^(3/2);

        dg5deXe = abs(N).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(N) + ...
            abs(N-1).*(5*e^2/2/(1-e^2)*b.^abs(N-1) + ...
            abs(N-2).*(e/sqrt(1-e^2)*b.^abs(N-2) + ...
            abs(N-3)/6.*b.^abs(N-3))))/sqrt(1-e^2) + ...
            abs(N-1).*abs(N)*e/(1-e^2)^(3/2).*b.^abs(N-2) +...
            5*abs(N)*e^2/(1-e^2)^2.*b.^abs(N-1) +...
            15*e^3/2/(1-e^2)^(5/2)*b.^abs(N);

        Jn = besselj(n,-k*e);
        % Jn/e nonsingular
        if n~=0
            Jn_e = -k/2/n*(besselj(n+1,-k*e)+besselj(n-1,-k*e));
        else
            Jn_e = Jn/e;
        end
        % dJnde/e nonsingular
        if abs(n)~=1
            dJnde_e = k^2/4/(n^2-1)*(2*Jn + (n+1)*besselj(n-2,-k*e) - (n-1)*besselj(n+2,-k*e));
        else
            dJnde_e = -0.5*k/e*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
        end
        % Jn/e^2 nonsingular
        if abs(n)>1
            Jn_e2 = k^2/4/n/(n^2-1)*(2*n*Jn + (n+1)*besselj(n-2,-k*e) + (n-1)*besselj(n+2,-k*e));
        else
            Jn_e2 = Jn/e^2;
        end

        dAk = Jn*a2.'*g5;
        dAk_e = Jn_e*a2.'*g5;
        dAkde_e = dJnde_e*a2.'*g5 + Jn_e*da2de.'*g5 + Jn_e2*a2.'*dg5deXe;

        dBk = Jn*b1.'*g5;
        dBk_e = Jn_e*b1.'*g5;
        dBkde_e = dJnde_e*b1.'*g5 + Jn_e*db1de.'*g5 + Jn_e2*b1.'*dg5deXe;

        dCk = Jn*g3;
        dCk_e = Jn_e*g3;
        dCkde_e = dJnde_e*g3  + Jn_e2*dg3deXe;

        % negative n
        n = -n;
        N = m + n + k + 1;
        Nc = n + k + 1;

        g3 = abs(Nc).*b.^abs(Nc-1) + e/sqrt(1-e^2)*b.^abs(Nc);

        g5 = abs(N-2).*abs(N-1).*abs(N)/6.*b.^abs(N-3) + ...
            e*abs(N-1).*abs(N)/sqrt(1-e^2).*b.^abs(N-2) + ...
            5*e^2*abs(N)/2/(1-e^2).*b.^abs(N-1) + ...
            5*e^3/2/(1-e^2)^(3/2).*b.^abs(N);

        dg3deXe = abs(Nc).*(abs(Nc-1).*b.^(abs(Nc-1)) + ...
            e/sqrt(1-e^2)*b.^(abs(Nc)))/sqrt(1-e^2) + ...
            e*b.^abs(Nc)/(1-e^2)^(3/2);

        dg5deXe = abs(N).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(N) + ...
            abs(N-1).*(5*e^2/2/(1-e^2)*b.^abs(N-1) + ...
            abs(N-2).*(e/sqrt(1-e^2)*b.^abs(N-2) + ...
            abs(N-3)/6.*b.^abs(N-3))))/sqrt(1-e^2) + ...
            abs(N-1).*abs(N)*e/(1-e^2)^(3/2).*b.^abs(N-2) +...
            5*abs(N)*e^2/(1-e^2)^2.*b.^abs(N-1) +...
            15*e^3/2/(1-e^2)^(5/2)*b.^abs(N);

        Jn = besselj(n,-k*e);
        % Jn/e nonsingular
        if n~=0
            Jn_e = -k/2/n*(besselj(n+1,-k*e)+besselj(n-1,-k*e));
        else
            Jn_e = Jn/e;
        end
        % dJnde/e nonsingular
        if abs(n)~=1
            dJnde_e = k^2/4/(n^2-1)*(2*Jn + (n+1)*besselj(n-2,-k*e) - (n-1)*besselj(n+2,-k*e));
        else
            dJnde_e = -0.5*k/e*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
        end
        % Jn/e^2 nonsingular
        if abs(n)>1
            Jn_e2 = k^2/4/n/(n^2-1)*(2*n*Jn + (n+1)*besselj(n-2,-k*e) + (n-1)*besselj(n+2,-k*e));
        else
            Jn_e2 = Jn/e^2;
        end

        dAk = dAk + Jn*a2.'*g5;
        dAk_e = dAk_e + Jn_e*a2.'*g5;
        dAkde_e = dAkde_e + dJnde_e*a2.'*g5 + Jn_e*da2de.'*g5 + Jn_e2*a2.'*dg5deXe;

        dBk = dBk + Jn*b1.'*g5;
        dBk_e = dBk_e + Jn_e*b1.'*g5;
        dBkde_e = dBkde_e + dJnde_e*b1.'*g5 + Jn_e*db1de.'*g5 + Jn_e2*b1.'*dg5deXe;

        dCk = dCk + Jn*g3;
        dCk_e = dCk_e + Jn_e*g3;
        dCkde_e = dCkde_e + dJnde_e*g3  + Jn_e2*dg3deXe;

        Ak = Ak + dAk;
        Ak_e = Ak_e + dAk_e;
        Akde_e = Akde_e + dAkde_e;

        Bk = Bk + dBk;
        Bk_e = Bk_e + dBk_e;
        Bkde_e = Bkde_e + dBkde_e;

        Ck = Ck + dCk;
        Ck_e = Ck_e + dCk_e;
        Ckde_e = Ckde_e + dCkde_e;

        n = abs(n); % return n to positive value
        n = n+1;
    end

    AkM(:,k) = 1/2/(1-e^2)^2*Ak;
    Ak_eM(:,k) = 1/2/(1-e^2)^2*Ak_e;
    Akde_eM(:,k) = 2/(1-e^2)^3*Ak + 1/2/(1-e^2)^2*Akde_e;

    BkM(:,k) = -1/2/(1-e^2)^1.5*Bk;
    Bk_eM(:,k) = -1/2/(1-e^2)^1.5*Bk_e;
    Bkde_eM(:,k) = -3/2/(1-e^2)^2.5*Bk - 1/2/(1-e^2)^1.5*Bkde_e;
    
    CkM(:,k) = 2/(1-e^2)*Ck;
    Ck_eM(:,k) = 2/(1-e^2)*Ck_e;
    Ckde_eM(:,k) = 4/(1-e^2)^2*Ck + 2/(1-e^2)*Ckde_e;
    
    k = k+1;
end


eta = sqrt(1-e^2);
% constant potential values
% R = -nMo^2*Re^2/2; % common factor
Rsma = -nMo*J2*Re^2/a;
Recc = -nMo*eta*J2*Re^2/2/a^2;
Rinc = -nMo*J2*Re^2*cos(i)/2/a^2/eta;
Rran = -nMo*J2*Re^2/2/a^2/eta;
Raop = -nMo*J2*Re^2/2/a^2;
Rman = Raop;

% Freq 0 elements without common factor
ran0 = 3*cos(i)/eta^3;
aop0 = -1.5*(5*cos(i)^2-1)/eta^4;
man0 = -1.5*(3*cos(i)^2-1)/eta^3;

% Calculate Spectrum of Elements
k = 1:kMax;
lpeSpec = nan(12,kMax);

lpeSpec(1:2,:) = Rsma*[s1*(BkM.*k);
    -(c1*AkM.*k+c0*CkM.*k)];
lpeSpec(3:4,:) = Recc*[eta*s1*(Bk_eM.*k) - (dc1do*Ak_eM+dc0do*Ck_eM);
    -eta*(c1*Ak_eM.*k+c0*Ck_eM.*k) - ds1do*Bk_eM];
lpeSpec(5:6,:) = Rinc*[dc1do_si*AkM + dc0do_si*CkM;
    ds1do_si*BkM];
lpeSpec(7:8,:) = Rran*[dc1di_si*AkM + dc0di_si*CkM;
    ds1di_si*BkM];
lpeSpec(9:10,:) = Raop*[eta*(c1*Akde_eM + c0*Ckde_eM) - cos(i)/eta*(dc1di_si*AkM + dc0di_si*CkM);
    eta*(s1*Bkde_eM) - cos(i)/eta*ds1di_si*BkM];
lpeSpec(11:12,:) = Rman*[-eta^2*(c1*Akde_eM + c0*Ckde_eM) + 6*(c1*AkM + c0*CkM);
    -eta^2*(s1*Bkde_eM) + 6*s1*BkM];
%% Zero Frequency Components
freq0 = [0; 0; 0; Rran*ran0; Raop*aop0; Rman*man0];
end