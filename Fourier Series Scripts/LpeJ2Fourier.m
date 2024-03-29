function [freq0, lpeSpec] = LpeJ2Fourier(oe,kMax,primary)

if nargin < 3
    primary = earth();
end

if size(oe,1) ~=6 && size(oe,2) ==6
    oe = oe.';
end
J2 = primary.J2;
Re = primary.Re;
mu = primary.mu;

% handle elements vector
sma = oe(1,:);
e = oe(2,:);
i = oe(3,:);
raan = oe(4,:);
aop = oe(5,:);
M = oe(6,:);

nT = size(oe,2);

nMo = sqrt(mu./sma.^3);
p = sma.*(1-e.^2);
b = (1-sqrt(1-e.^2))./e;

%% constant vectors
m2 = (0:4).';
m3 = (0:6).';
m4 = (0:8).';
m5 = (0:10).';

a2 = [ones(1,nT);-4*e;4*e.^2+2;-4*e;ones(1,nT)];
a3 = [ones(1,nT);-6*e;12*e.^2+3;-8*e.^3-12*e;12*e.^2+3;-6*e;ones(1,nT)];
a4 = [ones(1,nT);-8*e;24*e.^2+4;-(32*e.^3+24*e);(16*e.^4+48*e.^2+6);...
    -(32*e.^3+24*e);24*e.^2+4;-8*e;ones(1,nT)];
a5 = [ones(1,nT);-10*e;40*e.^2+5;-(80*e.^3+40*e);(80*e.^4+120*e.^2+10);...
    -(32*e.^5+160*e.^3+60*e);(80*e.^4+120*e.^2+10);-(80*e.^3+40*e);40*e.^2+5;...
    -10*e;ones(1,nT)];

b1 = [-ones(1,nT);2*e;zeros(1,nT);-2*e;ones(1,nT)];
b2 = [-ones(1,nT);4*e;-4*e.^2-1;zeros(1,nT);4*e.^2+1;-4*e;ones(1,nT)];
b3 = [-ones(1,nT);6*e;-12*e.^2-2;8*e.^3+6*e;zeros(1,nT);-8*e.^3-6*e;12*e.^2+2;...
    -6*e;ones(1,nT)];
b4 = [-ones(1,nT);8*e;-24*e.^2-3;32*e.^3+16*e;-16*e.^4-24*e.^2-2;zeros(1,nT);...
    16*e.^4+24*e.^2+2;-32*e.^3-16*e;24*e.^2+3;-8*e;ones(1,nT)];

da2de = [zeros(1,nT);-4*ones(1,nT);8*e;-4*ones(1,nT);zeros(1,nT)];
da3de = [zeros(1,nT);-6*ones(1,nT);24*e;(-24*e.^2-12);24*e;-6*ones(1,nT);zeros(1,nT)];
da4de = [zeros(1,nT);-8*ones(1,nT);48*e;(-96*e.^2-24);(64*e.^3+96*e);...
    (-96*e.^2-24);48*e;-8*ones(1,nT);zeros(1,nT)];
da5de = [zeros(1,nT);-10*ones(1,nT);80*e;(-240*e.^2-40);(320*e.^3+240*e);...
    (-160*e.^4-480*e.^2-60);(320*e.^3+240*e);(-240*e.^2-40);80*e;-10*ones(1,nT);zeros(1,nT)];

db1de = [zeros(1,nT);2*ones(1,nT);zeros(1,nT);-2*ones(1,nT);zeros(1,nT)];
db2de = [zeros(1,nT);4*ones(1,nT);-8*e;zeros(1,nT);8*e;-4*ones(1,nT);zeros(1,nT)];
db3de = [zeros(1,nT);6*ones(1,nT);-24*e;24*e.^2+6;zeros(1,nT);-24*e.^2-6;24*e;...
    -6*ones(1,nT);zeros(1,nT)];
db4de = [zeros(1,nT);8*ones(1,nT);-48*e;96*e.^2+16;-64*e.^3-48*e;zeros(1,nT);...
    64*e.^3+48*e;-96*e.^2-16;48*e;-8*ones(1,nT);zeros(1,nT)];

%% matrices by aop
C = [6*(3*sin(i).^2.*cos(aop).^2-1)./(1-e.^2).^2;
    (9*e.^2.*sin(i).^2.*cos(aop).^2+3*sin(i).^2.*(sin(aop).^2-cos(aop).^2)-3*e.^2)/2./(1-e.^2).^3.5;
    (3*e.^3.*sin(i).^2.*cos(aop).^2+9*e.*sin(i).^2.*(sin(aop).^2-cos(aop).^2)-e.^3)/4./(1-e.^2).^4;
    (9*e.^2.*sin(i).^2.*(sin(aop).^2-cos(aop).^2))/8./(1-e.^2).^4.5;
    (3*e.^3.*sin(i).^2.*(sin(aop).^2-cos(aop).^2))/16./(1-e.^2).^5];

dCdi_si = [36*cos(aop).^2.*cos(i)./(1-e.^2).^2; % all pre-divided by sin(i)
    3*cos(i).*(3*e.^2.*cos(aop).^2-2*cos(aop).^2+1)./(1-e.^2).^3.5;
    3*e.*cos(i).*(3+(e.^2-6).*cos(aop).^2)/2./(1-e.^2).^4;
    -9*e.^2.*cos(i).*(2*cos(aop).^2-1)/4./(1-e.^2).^4.5;
    -3*e.^3.*cos(i).*(2*cos(aop).^2-1)/8./(1-e.^2).^5];

dCdo = [-36*sin(i).^2.*sin(aop).*cos(aop)./(1-e.^2).^2;
    3*(-3*e.^2+2).*sin(i).^2.*sin(aop).*cos(aop)./(1-e.^2).^3.5;
    -3*e.*(e.^2-6).*sin(i).^2.*sin(aop).*cos(aop)/2./(1-e.^2).^4;
    9*e.^2.*sin(i).^2.*sin(aop).*cos(aop)/2./(1-e.^2).^4.5;
    3*e.^3.*sin(i).^2.*sin(aop).*cos(aop)/4./(1-e.^2).^5];

dCdo_si =[-36*sin(i).*sin(aop).*cos(aop)./(1-e.^2).^2;
    3*(-3*e.^2+2).*sin(i).*sin(aop).*cos(aop)./(1-e.^2).^3.5;
    -3*e.*(e.^2-6).*sin(i).*sin(aop).*cos(aop)/2./(1-e.^2).^4;
    9*e.^2.*sin(i).*sin(aop).*cos(aop)/2./(1-e.^2).^4.5;
    3*e.^3.*sin(i).*sin(aop).*cos(aop)/4./(1-e.^2).^5];

dCde = [24*e.*(3*sin(i).^2.*cos(aop).^2-1)./(1-e.^2).^3;
    3*e.*(15*e.^2.*sin(i).^2.*cos(aop).^2 +7*sin(i).^2.*sin(aop).^2 ...
    -sin(i).^2.*cos(aop).^2 -5*e.^2 -2)/2./(1-e.^2).^4.5;
    ((15*e.^4-54*e.^2-9).*sin(i).^2.*cos(aop).^2 + ...
    (63*e.^2+9).*sin(i).^2.*sin(aop).^2-5*e.^4-3*e.^2)/4./(1-e.^2).^5;
    -9*e.*sin(i).^2.*(cos(aop).^2-sin(aop).^2).*(7*e.^2+2)/8./(1-e.^2).^5.5;
    -3*e.^2.*sin(i).^2.*(cos(aop).^2-sin(aop).^2).*(7*e.^2+3)/16./(1-e.^2).^6];

S = -[6*sin(i).^2.*sin(aop).*cos(aop)/2./(1-e.^2).^3;
    18*e.*sin(i).^2.*sin(aop).*cos(aop)/4./(1-e.^2).^3.5;
    18*e.^2.*sin(i).^2.*sin(aop).*cos(aop)/8./(1-e.^2).^4;
    6*e.^3.*sin(i).^2.*sin(aop).*cos(aop)/16./(1-e.^2).^4.5];

dSdi_si = -[6*cos(i).*sin(aop).*cos(aop)./(1-e.^2).^3;
    9*e.*cos(i).*sin(aop).*cos(aop)./(1-e.^2).^3.5;
    9*e.^2.*cos(i).*sin(aop).*cos(aop)/2./(1-e.^2).^4;
    3*e.^3.*cos(i).*sin(aop).*cos(aop)/4./(1-e.^2).^4.5];

dSdo = -[6*sin(i).^2.*(2*cos(aop).^2-1)/2./(1-e.^2).^3;
    18*e.*sin(i).^2.*(2*cos(aop).^2-1)/4./(1-e.^2).^3.5;
    18*e.^2.*sin(i).^2.*(2*cos(aop).^2-1)/8./(1-e.^2).^4;
    6*e.^3.*sin(i).^2.*(2*cos(aop).^2-1)/16./(1-e.^2).^4.5];

dSdo_si = -[6*sin(i).*(2*cos(aop).^2-1)/2./(1-e.^2).^3;
    18*e.*sin(i).*(2*cos(aop).^2-1)/4./(1-e.^2).^3.5;
    18*e.^2.*sin(i).*(2*cos(aop).^2-1)/8./(1-e.^2).^4;
    6*e.^3.*sin(i).*(2*cos(aop).^2-1)/16./(1-e.^2).^4.5];

dSde = [-18*e.*sin(i).^2.*sin(aop).*cos(aop)./(1-e.^2).^4;
    -9*sin(i).^2.*sin(aop).*cos(aop).*(6*e.^2+1)/2./(1-e.^2).^4.5;
    -9*e.*sin(i).^2.*sin(aop).*cos(aop).*(3*e.^2+1)/2./(1-e.^2).^5;
    -9*e.^2.*sin(i).^2.*sin(aop).*cos(aop).*(2*e.^2+1)/8./(1-e.^2).^5.5];

%% Second Order - O(J2^2), not necessary

% d2Cdo2 = [-36*sin(i)^2*(2*cos(aop).^2-1)/(1-e^2)^2;
%     -3*(3*e^2-2)*sin(i)^2*(2*cos(aop).^2-1)/(1-e^2)^3.5;
%     -3*(e^2-6)*e*sin(i)^2*(2*cos(aop).^2-1)/2/(1-e^2)^4;
%     9*e^2*sin(i)^2*(2*cos(aop).^2-1)/2/(1-e^2)^4.5;
%     3*e^3*sin(i)^2*(2*cos(aop).^2-1)/4/(1-e^2)^5];
%
% d2Cdo2_si = [-36*sin(i)*(2*cos(aop).^2-1)/(1-e^2)^2;
%     -3*(3*e^2-2)*sin(i)*(2*cos(aop).^2-1)/(1-e^2)^3.5;
%     -3*(e^2-6)*e*sin(i)*(2*cos(aop).^2-1)/2/(1-e^2)^4;
%     9*e^2*sin(i)*(2*cos(aop).^2-1)/2/(1-e^2)^4.5;
%     3*e^3*sin(i)*(2*cos(aop).^2-1)/4/(1-e^2)^5];
%
% d2Cdido_si = [-72*cos(i)*sin(aop).*cos(aop)/(1-e^2)^2;
%     (12-18*e^2)*cos(i)*sin(aop).*cos(aop)/(1-e^2)^3.5;
%     3*e*(6-e^2)*cos(i)*sin(aop).*cos(aop)/(1-e^2)^4;
%     9*e^2*cos(i)*sin(aop).*cos(aop)/(1-e^2)^4.5;
%     3*e^3*cos(i)*sin(aop).*cos(aop)/2/(1-e^2)^5];
%
% d2Cdedo = [-144*e*sin(i)^2*sin(aop).*cos(aop)/(1-e^2)^3;
%     -(45*e^2-24)*e*sin(i)^2*sin(aop).*cos(aop)/(1-e^2)^4.5;
%     -3*(5*e^4-39*e^2-6)*sin(i)^2*sin(aop).*cos(aop)/2/(1-e^2)^5;
%     (63*e^2+18)*e*sin(i)^2*sin(aop).*cos(aop)/2/(1-e^2)^5.5;
%     -3*e^2*(7*e^2+3)*sin(i)^2*sin(aop).*cos(aop)/4/(1-e^2)^6];
%
% d2Sdo2 = [12*sin(i)^2*sin(aop).*cos(aop)/(1-e^2)^3;
%     18*e*sin(i)^2*sin(aop).*cos(aop)/(1-e^2)^3.5;
%     9*e^2*sin(i)^2*sin(aop).*cos(aop)/(1-e^2)^4;
%     3*e^3*sin(i)^2*sin(aop).*cos(aop)/2/(1-e^2)^4.5];
%
% d2Sdo2_si = [12*sin(i)*sin(aop).*cos(aop)/(1-e^2)^3;
%     18*e*sin(i)*sin(aop).*cos(aop)/(1-e^2)^3.5;
%     9*e^2*sin(i)*sin(aop).*cos(aop)/(1-e^2)^4;
%     3*e^3*sin(i)*sin(aop).*cos(aop)/2/(1-e^2)^4.5];
%
% d2Sdido_si = [-6*cos(i)*(2*cos(aop).^2-1)/(1-e^2)^3;
%     -9*e*cos(i)*(2*cos(aop).^2-1)/(1-e^2)^3.5;
%     -9*e^2*cos(i)*(2*cos(aop).^2-1)/2/(1-e^2)^4;
%     -3*e^3*cos(i)*(2*cos(aop).^2-1)/4/(1-e^2)^4.5];
%
% d2Sdedo = [-18*e*sin(i)^2*(2*cos(aop).^2-1)/(1-e^2)^4;
%     -9*(6*e^2+1)*sin(i)^2*(2*cos(aop).^2-1)/2/(1-e^2)^4.5;
%     -9*e*(3*e^2+1)*sin(i)^2*(2*cos(aop).^2-1)/2/(1-e^2)^5;
%     -9*e^2*(2*e^2+1)*sin(i)^2*(2*cos(aop).^2-1)/8/(1-e^2)^5.5];
%% Loop for A,B
AkM = nan(5,nT,kMax);
Ak_eM = nan(5,nT,kMax);
Akde_eM = nan(5,nT,kMax);

BkM = nan(4,nT,kMax);
Bk_eM = nan(4,nT,kMax);
Bkde_eM = nan(4,nT,kMax);

k = 1;

while k <= kMax
    n = 0;

    g2 = b.^abs(m2+n+k-2);
    g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e./sqrt(1-e.^2).*b.^abs(m3+n+k-2);
    g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
        3*e.*abs(m4+n+k-2)/2./sqrt(1-e.^2).*b.^abs(m4+n+k-3) + ...
        3*e.^2/2./(1-e.^2).*b.^abs(m4+n+k-2);
    g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
        e.*abs(m5+n+k-3).*abs(m5+n+k-2)./sqrt(1-e.^2).*b.^abs(m5+n+k-4) + ...
        5*e.^2.*abs(m5+n+k-2)/2./(1-e.^2).*b.^abs(m5+n+k-3) + ...
        5*e.^3/2./(1-e.^2).^(3/2).*b.^abs(m5+n+k-2);

    dg2deXe = abs(m2+n+k-2).*b.^(abs(m2+n+k-2))./sqrt(1-e.^2);
    dg3deXe = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^(abs(m3+n+k-3)) + ...
        e./sqrt(1-e.^2).*b.^(abs(m3+n+k-2)))./sqrt(1-e.^2) + ...
        e.*b.^abs(m3+n+k-2)./(1-e.^2).^(3/2);
    dg4deXe = abs(m4+n+k-2).*(3*e.^2/2./(1-e.^2).*b.^abs(m4+n+k-2) +...
        abs(m4+n+k-3).*(3*e/2./sqrt(1-e.^2).*b.^abs(m4+n+k-3) +...
        abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))./sqrt(1-e.^2) +...
        3/2*e.*abs(m4+n+k-2)./(1-e.^2).^(3/2).*b.^abs(m4+n+k-3) +...
        3*e.^2./(1-e.^2).^2.*b.^abs(m4+n+k-2);
    dg5deXe = abs(m5+n+k-2).*(5*e.^3/2./(1-e.^2).^(3/2).*b.^abs(m5+n+k-2) + ...
        abs(m5+n+k-3).*(5*e.^2/2./(1-e.^2).*b.^abs(m5+n+k-3) + ...
        abs(m5+n+k-4).*(e./sqrt(1-e.^2).*b.^abs(m5+n+k-4) + ...
        abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))./sqrt(1-e.^2) + ...
        abs(m5+n+k-3).*abs(m5+n+k-2).*e./(1-e.^2).^(3/2).*b.^abs(m5+n+k-4) +...
        5*abs(m5+n+k-2).*e.^2./(1-e.^2).^2.*b.^abs(m5+n+k-3) +...
        15*e.^3/2./(1-e.^2).^(5/2).*b.^abs(m5+n+k-2);

    Jk = besselj(k,k*e);
    % Jk/e nonsingular
    Jk_e = 0.5*(besselj(k+1,k*e) + besselj(k-1,k*e));
    % dJkde/e nonsingular
    if k~=1
        % use expression with no /e
        dJkde_e = k^2/4/(k^2-1)*(2*Jk + (k+1)*besselj(k-2,k*e) - (k-1)*besselj(k+2,k*e));
    else
        % elimination of /e not possible
        dJkde_e = 0.5*k./e.*(besselj(k-1,k*e) - besselj(k+1,k*e));
    end

    Jn = besselj(n,-k*e);
    % Jn/e nonsingular
    if n~=0
        Jn_e = -k/2/n*(besselj(n+1,-k*e)+besselj(n-1,-k*e));
    else
        Jn_e = Jn./e;
    end
    % dJnde/e nonsingular
    if abs(n)~=1
        dJnde_e = k^2/4/(n^2-1)*(2*Jn + (n+1)*besselj(n-2,-k*e) - (n-1)*besselj(n+2,-k*e));
    else
        dJnde_e = -0.5*k./e.*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
    end
    % Jn/e^2 nonsingular
    if abs(n)>1
        Jn_e2 = k^2/4/n/(n^2-1)*(2*n*Jn + (n+1)*besselj(n-2,-k*e) + (n-1)*besselj(n+2,-k*e));
    else
        Jn_e2 = Jn./e.^2;
    end

    Ak = [Jk; Jn.*[dot(a2,g2); dot(a3,g3); dot(a4,g4); dot(a5,g5)]];
    Ak_e = [Jk_e; Jn_e.*[dot(a2,g2); dot(a3,g3); dot(a4,g4); dot(a5,g5)]];
    Akde_e = [dJkde_e; dJnde_e.*[dot(a2,g2); dot(a3,g3); dot(a4,g4); dot(a5,g5)] +...
        Jn_e.*[dot(da2de,g2); dot(da3de,g3); dot(da4de,g4); dot(da5de,g5)] +...
        Jn_e2.*[dot(a2,dg2deXe); dot(a3,dg3deXe); dot(a4,dg4deXe); dot(a5,dg5deXe)]];

    Bk = Jn.*[dot(b1,g2); dot(b2,g3); dot(b3,g4); dot(b4,g5)];
    Bk_e = Jn_e.*[dot(b1,g2); dot(b2,g3); dot(b3,g4); dot(b4,g5)];
    Bkde_e = dJnde_e.*[dot(b1,g2); dot(b2,g3); dot(b3,g4); dot(b4,g5)] +...
        Jn_e.*[dot(db1de,g2); dot(db2de,g3); dot(db3de,g4); dot(db4de,g5)] +...
        Jn_e2.*[dot(b1,dg2deXe); dot(b2,dg3deXe); dot(b3,dg4deXe); dot(b4,dg5deXe)];

    n = 1;
    while n <= k + 5
        % positive n
        g2 = b.^abs(m2+n+k-2);
        g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e./sqrt(1-e.^2).*b.^abs(m3+n+k-2);
        g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
            3*e.*abs(m4+n+k-2)/2./sqrt(1-e.^2).*b.^abs(m4+n+k-3) + ...
            3*e.^2/2./(1-e.^2).*b.^abs(m4+n+k-2);
        g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
            e.*abs(m5+n+k-3).*abs(m5+n+k-2)./sqrt(1-e.^2).*b.^abs(m5+n+k-4) + ...
            5*e.^2.*abs(m5+n+k-2)/2./(1-e.^2).*b.^abs(m5+n+k-3) + ...
            5*e.^3/2./(1-e.^2).^(3/2).*b.^abs(m5+n+k-2);

        dg2deXe = abs(m2+n+k-2).*b.^(abs(m2+n+k-2))./sqrt(1-e.^2);
        dg3deXe = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^(abs(m3+n+k-3)) + ...
            e./sqrt(1-e.^2).*b.^(abs(m3+n+k-2)))./sqrt(1-e.^2) + ...
            e.*b.^abs(m3+n+k-2)./(1-e.^2).^(3/2);
        dg4deXe = abs(m4+n+k-2).*(3*e.^2/2./(1-e.^2).*b.^abs(m4+n+k-2) +...
            abs(m4+n+k-3).*(3*e/2./sqrt(1-e.^2).*b.^abs(m4+n+k-3) +...
            abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))./sqrt(1-e.^2) +...
            3/2*e.*abs(m4+n+k-2)./(1-e.^2).^(3/2).*b.^abs(m4+n+k-3) +...
            3*e.^2./(1-e.^2).^2.*b.^abs(m4+n+k-2);
        dg5deXe = abs(m5+n+k-2).*(5*e.^3/2./(1-e.^2).^(3/2).*b.^abs(m5+n+k-2) + ...
            abs(m5+n+k-3).*(5*e.^2/2./(1-e.^2).*b.^abs(m5+n+k-3) + ...
            abs(m5+n+k-4).*(e./sqrt(1-e.^2).*b.^abs(m5+n+k-4) + ...
            abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))./sqrt(1-e.^2) + ...
            abs(m5+n+k-3).*abs(m5+n+k-2).*e./(1-e.^2).^(3/2).*b.^abs(m5+n+k-4) +...
            5*abs(m5+n+k-2).*e.^2./(1-e.^2).^2.*b.^abs(m5+n+k-3) +...
            15*e.^3/2./(1-e.^2).^(5/2).*b.^abs(m5+n+k-2);

        Jn = besselj(n,-k*e);
        % Jn/e nonsingular
        if n~=0
            Jn_e = -k/2/n*(besselj(n+1,-k*e)+besselj(n-1,-k*e));
        else
            Jn_e = Jn./e;
        end
        % dJnde/e nonsingular
        if abs(n)~=1
            dJnde_e = k^2/4/(n^2-1)*(2*Jn + (n+1)*besselj(n-2,-k*e) - (n-1)*besselj(n+2,-k*e));
        else
            dJnde_e = -0.5*k./e.*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
        end
        % Jn/e^2 nonsingular
        if abs(n)>1
            Jn_e2 = k^2/4/n/(n^2-1)*(2*n*Jn + (n+1)*besselj(n-2,-k*e) + (n-1)*besselj(n+2,-k*e));
        else
            Jn_e2 = Jn./e.^2;
        end

        dAk = [zeros(1,nT); Jn.*[dot(a2,g2); dot(a3,g3); dot(a4,g4); dot(a5,g5)]];
        dAk_e = [zeros(1,nT); Jn_e.*[dot(a2,g2); dot(a3,g3); dot(a4,g4); dot(a5,g5)]];
        dAkde_e = [zeros(1,nT); dJnde_e.*[dot(a2,g2); dot(a3,g3); dot(a4,g4); dot(a5,g5)] +...
            Jn_e.*[dot(da2de,g2); dot(da3de,g3); dot(da4de,g4); dot(da5de,g5)] +...
            Jn_e2.*[dot(a2,dg2deXe); dot(a3,dg3deXe); dot(a4,dg4deXe); dot(a5,dg5deXe)]];

        dBk = Jn.*[dot(b1,g2); dot(b2,g3); dot(b3,g4); dot(b4,g5)];
        dBk_e = Jn_e.*[dot(b1,g2); dot(b2,g3); dot(b3,g4); dot(b4,g5)];
        dBkde_e = dJnde_e.*[dot(b1,g2); dot(b2,g3); dot(b3,g4); dot(b4,g5)] +...
            Jn_e.*[dot(db1de,g2); dot(db2de,g3); dot(db3de,g4); dot(db4de,g5)] +...
            Jn_e2.*[dot(b1,dg2deXe); dot(b2,dg3deXe); dot(b3,dg4deXe); dot(b4,dg5deXe)];


        % negative n
        n = -n;
        g2 = b.^abs(m2+n+k-2);
        g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e./sqrt(1-e.^2).*b.^abs(m3+n+k-2);
        g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
            3*e.*abs(m4+n+k-2)/2./sqrt(1-e.^2).*b.^abs(m4+n+k-3) + ...
            3*e.^2/2./(1-e.^2).*b.^abs(m4+n+k-2);
        g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
            e.*abs(m5+n+k-3).*abs(m5+n+k-2)./sqrt(1-e.^2).*b.^abs(m5+n+k-4) + ...
            5*e.^2.*abs(m5+n+k-2)/2./(1-e.^2).*b.^abs(m5+n+k-3) + ...
            5*e.^3/2./(1-e.^2).^(3/2).*b.^abs(m5+n+k-2);

        dg2deXe = abs(m2+n+k-2).*b.^(abs(m2+n+k-2))./sqrt(1-e.^2);
        dg3deXe = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^(abs(m3+n+k-3)) + ...
            e./sqrt(1-e.^2).*b.^(abs(m3+n+k-2)))./sqrt(1-e.^2) + ...
            e.*b.^abs(m3+n+k-2)./(1-e.^2).^(3/2);
        dg4deXe = abs(m4+n+k-2).*(3*e.^2/2./(1-e.^2).*b.^abs(m4+n+k-2) +...
            abs(m4+n+k-3).*(3*e/2./sqrt(1-e.^2).*b.^abs(m4+n+k-3) +...
            abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))./sqrt(1-e.^2) +...
            3/2*e.*abs(m4+n+k-2)./(1-e.^2).^(3/2).*b.^abs(m4+n+k-3) +...
            3*e.^2./(1-e.^2).^2.*b.^abs(m4+n+k-2);
        dg5deXe = abs(m5+n+k-2).*(5*e.^3/2./(1-e.^2).^(3/2).*b.^abs(m5+n+k-2) + ...
            abs(m5+n+k-3).*(5*e.^2/2./(1-e.^2).*b.^abs(m5+n+k-3) + ...
            abs(m5+n+k-4).*(e./sqrt(1-e.^2).*b.^abs(m5+n+k-4) + ...
            abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))./sqrt(1-e.^2) + ...
            abs(m5+n+k-3).*abs(m5+n+k-2).*e./(1-e.^2).^(3/2).*b.^abs(m5+n+k-4) +...
            5*abs(m5+n+k-2).*e.^2./(1-e.^2).^2.*b.^abs(m5+n+k-3) +...
            15*e.^3/2./(1-e.^2).^(5/2).*b.^abs(m5+n+k-2);

        Jn = besselj(n,-k*e);
        % Jn/e nonsingular
        if n~=0
            Jn_e = -k/2/n*(besselj(n+1,-k*e)+besselj(n-1,-k*e));
        else
            Jn_e = Jn./e;
        end
        % dJnde/e nonsingular
        if abs(n)~=1
            dJnde_e = k^2/4/(n^2-1)*(2*Jn + (n+1)*besselj(n-2,-k*e) - (n-1)*besselj(n+2,-k*e));
        else
            dJnde_e = -0.5*k./e.*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
        end
        % Jn/e^2 nonsingular
        if abs(n)>1
            Jn_e2 = k^2/4/n/(n^2-1)*(2*n*Jn + (n+1)*besselj(n-2,-k*e) + (n-1)*besselj(n+2,-k*e));
        else
            Jn_e2 = Jn./e.^2;
        end

        dAk = dAk + [zeros(1,nT); Jn.*[dot(a2,g2); dot(a3,g3); dot(a4,g4); dot(a5,g5)]];
        dAk_e = dAk_e + [zeros(1,nT); Jn_e.*[dot(a2,g2); dot(a3,g3); dot(a4,g4); dot(a5,g5)]];
        dAkde_e = dAkde_e + [zeros(1,nT); dJnde_e.*[dot(a2,g2); dot(a3,g3); dot(a4,g4); dot(a5,g5)] +...
            Jn_e.*[dot(da2de,g2); dot(da3de,g3); dot(da4de,g4); dot(da5de,g5)] +...
            Jn_e2.*[dot(a2,dg2deXe); dot(a3,dg3deXe); dot(a4,dg4deXe); dot(a5,dg5deXe)]];

        dBk = dBk + Jn.*[dot(b1,g2); dot(b2,g3); dot(b3,g4); dot(b4,g5)];
        dBk_e = dBk_e + Jn_e.*[dot(b1,g2); dot(b2,g3); dot(b3,g4); dot(b4,g5)];
        dBkde_e = dBkde_e + dJnde_e.*[dot(b1,g2); dot(b2,g3); dot(b3,g4); dot(b4,g5)] +...
            Jn_e.*[dot(db1de,g2); dot(db2de,g3); dot(db3de,g4); dot(db4de,g5)] +...
            Jn_e2.*[dot(b1,dg2deXe); dot(b2,dg3deXe); dot(b3,dg4deXe); dot(b4,dg5deXe)];

        Ak = Ak + dAk;
        Ak_e = Ak_e + dAk_e;
        Akde_e = Akde_e + dAkde_e;

        Bk = Bk + dBk;
        Bk_e = Bk_e + dBk_e;
        Bkde_e = Bkde_e + dBkde_e;

        n = abs(n); % return n to positive value
        n = n+1;
    end

    AkM(:,:,k) = Ak;
    Ak_eM(:,:,k) = Ak_e;
    Akde_eM(:,:,k) = Akde_e;

    BkM(:,:,k) = Bk;
    Bk_eM(:,:,k) = Bk_e;
    Bkde_eM(:,:,k) = Bkde_e;

    k = k+1;
end
%% Define Final constants
eta = sqrt(1-e.^2);
% constant potential values
% R = -nMo^2*Re^2/2; % common factor
Rsma = -nMo*J2*Re^2./sma;
Recc = -nMo.*eta*J2*Re^2/2./sma.^2;
Rinc = -nMo*J2*Re^2.*cos(i)/2./sma.^2./eta;
Rran = -nMo*J2*Re^2/2./sma.^2./eta;
Raop = -nMo*J2*Re^2/2./sma.^2;
Rman = Raop;

% Freq 0 elements without common factor
ran0 = 3*cos(i)./eta.^3;
aop0 = -1.5*(5*cos(i).^2-1)./eta.^4;
man0 = -1.5*(3*cos(i).^2-1)./eta.^3;

% Calculate Spectrum of Elements
dSmaSpec = nan(2,nT,kMax);
dEccSpec = nan(2,nT,kMax);
dIncSpec = nan(2,nT,kMax);
dRanSpec = nan(2,nT,kMax);
dAopSpec = nan(2,nT,kMax);
dManSpec = nan(2,nT,kMax);

for k = 1:kMax

dSmaSpec(:,:,k) = Rsma.*[dot(S,(squeeze(BkM(:,:,k)))*k);...
    -dot(C,(squeeze(AkM(:,:,k)))*k)];

dEccSpec(:,:,k) = Recc.*[eta.*dot(S,squeeze(Bk_eM(:,:,k)))*k - dot(dCdo,squeeze(Ak_eM(:,:,k)));...
    -eta.*dot(C,squeeze(Ak_eM(:,:,k)))*k - dot(dSdo,squeeze(Bk_eM(:,:,k)))];

dIncSpec(:,:,k) = Rinc.*[dot(dCdo_si,squeeze(AkM(:,:,k)));...
    dot(dSdo_si,squeeze(BkM(:,:,k)))];

dRanSpec(:,:,k) = Rran.*[dot(dCdi_si,squeeze(AkM(:,:,k)));...
    dot(dSdi_si,squeeze(BkM(:,:,k)))];

dAopSpec(:,:,k) = Raop.*[eta.*(dot(dCde,squeeze(Ak_eM(:,:,k))) + dot(C,squeeze(Akde_eM(:,:,k))))...
    - cos(i)./eta.*dot(dCdi_si,squeeze(AkM(:,:,k)));...
    eta.*(dot(dSde,squeeze(Bk_eM(:,:,k))) + dot(S,squeeze(Bkde_eM(:,:,k))))...
    - cos(i)./eta.*dot(dSdi_si,squeeze(BkM(:,:,k)))];

dManSpec(:,:,k) = Rman.*[-eta.^2.*(dot(dCde,squeeze(Ak_eM(:,:,k))) + dot(C,squeeze(Akde_eM(:,:,k))))...
    + 6*dot(C,squeeze(AkM(:,:,k)));...
    -eta.^2.*(dot(dSde,squeeze(Bk_eM(:,:,k))) + dot(S,squeeze(Bkde_eM(:,:,k))))...
    + 6*dot(S,squeeze(BkM(:,:,k)))];

% dEccSpec = Recc*[eta*S.'*(Bk_eM.*k) - dCdo.'*Ak_eM,...
%     -eta*C.'*(Ak_eM.*k) - dSdo.'*Bk_eM];
% 
% dIncSpec = Rinc*[dCdo_si.'*AkM,...
%     dSdo_si.'*BkM];
% 
% dRanSpec = Rran*[dCdi_si.'*AkM,...
%     dSdi_si.'*BkM];
%
% dAopSpec = Raop*[eta*(dCde.'*Ak_eM + C.'*Akde_eM) - cos(i)/eta*dCdi_si.'*AkM,...
%     eta*(dSde.'*Bk_eM + S.'*Bkde_eM) - cos(i)/eta*dSdi_si.'*BkM];
% 
% dManSpec = Rman*[-eta^2*(dCde.'*Ak_eM + C.'*Akde_eM) + 6*C.'*AkM,...
%     -eta^2*(dSde.'*Bk_eM + S.'*Bkde_eM) + 6*S.'*BkM];

end
lpeSpec = [dSmaSpec;
    dEccSpec;
    dIncSpec;
    dRanSpec;
    dAopSpec;
    dManSpec];

%% Zero Frequency Components
freq0 = [zeros(3,nT); Rran.*ran0; Raop.*aop0; Rman.*man0];

