clear
% tol = 1e-14;
tol = 20;
primary = earth();
%% Orbital Elements
oe = [10000, 0.3, 20*pi/180, 0, 70*pi/180, 0];
a = oe(1);
e = oe(2);
i = oe(3);
raan = oe(4);
aop = oe(5);
M = 0:0.01:2*pi;
f = pi/180*me2ta(180/pi*M,e,1e-14);
%% Calculate Potentials

% [J2F,J2freq] = J2PotFourier(oe,M,tol);
%
% J2T = -primary.mu*primary.J2*primary.Re^2/2/a^3/(1-e^2)^3*(1+e*cos(f)).^3.*...
%     (3*sin(aop+f).^2*sin(i)^2-1);

dJ2diT = -3*primary.mu*primary.J2*primary.Re^2/a^3/(1-e^2)^3*(1+e*cos(f)).^3.*...
    (sin(f)*cos(aop)+cos(f)*sin(aop)).^2*sin(i)*cos(i);

dJ2doT = -3*primary.mu*primary.J2*primary.Re^2/a^3/(1-e^2)^3*(1+e*cos(f)).^3.*...
    sin(f+aop).*cos(f+aop)*sin(i)^2;

dJ2dlT = 3*primary.mu*primary.J2*primary.Re^2/2/a^3/(1-e^2)^(9/2)*(1+e*cos(f)).^4.*...
    (sin(i)^2*sin(f+aop).*(3*e*sin(f+aop).*sin(f) - 2*e*cos(f+aop).*cos(f)-2*cos(f+aop))-...
    e*sin(f));

dJ2daT = 3*primary.mu*primary.J2*primary.Re^2/2/a^4/(1-e^2)^3*(1+e*cos(f)).^3.*...
    (3*sin(aop+f).^2*sin(i)^2-1);

% eccentiricity true derivative
C = [9*e*sin(i)^2*cos(aop)^2-3*e;
    9*e^2*sin(i)^2*cos(aop)^2+3*sin(i)^2*(sin(aop)^2-cos(aop)^2)-3*e^2;
    3*e^3*sin(i)^2*cos(aop)^2+9*e*sin(i)^2*(sin(aop)^2-cos(aop)^2)-e^3;
    9*e^2*sin(i)^2*(sin(aop)^2-cos(aop)^2);
    3*e^3*sin(i)^2*(sin(aop)^2-cos(aop)^2)];

S = [6*sin(i)^2*sin(aop)*cos(aop);
    18*e*sin(i)^2*sin(aop)*cos(aop);
    18*e^2*sin(i)^2*sin(aop)*cos(aop);
    6*e^3*sin(i)^2*sin(aop)*cos(aop)];

dCDende = [3*(3*sin(i)^2*cos(aop)^2-1)*(5*e^2+1)/(1-e^2)^4;
           6*e*(6*e^2*sin(i)^2*cos(aop)^2+3*sin(i)^2*sin(aop)^2-2*e^2-1)/(1-e^2)^4;
           3*((3*e^4-12*e^2-3)*sin(i)^2*cos(aop)^2 + (15*e^2+3)*sin(i)^2*sin(aop)^2 -e^4-e^2)/(1-e^2)^4;
           18*e*sin(i)^2*(sin(aop)^2-cos(aop)^2)*(2*e^2+1)/(1-e^2)^4;
           9*e^2*sin(i)^2*(sin(aop)^2-cos(aop)^2)*(e^2+1)/(1-e^2)^4];
            
dSDende = [36*e*sin(i)^2*sin(aop)*cos(aop)/(1-e^2)^4;
    18*sin(i)^2*sin(aop)*cos(aop)*(5*e^2+1)/(1-e^2)^4;
    36*e*sin(i)^2*sin(aop)*cos(aop)*(2*e^2+1)/(1-e^2)^4;
    18*e^2*sin(i)^2*sin(aop)*cos(aop)*(e^2+1)/(1-e^2)^4];

dfde = (sin(f).*(e*cos(f)+2))/(1-e^2);

% ***** dRde separated by frequency, for debugging only *****
% dJ2deT = -primary.mu*primary.J2*primary.Re^2/2/a^3*(...
%     (6*e*(3*sin(i)^2*cos(aop)^2-1)/(1-e^2)^4 + ...
%     dCDende(1)*cos(f) -...
%     C(1)*sin(f).*dfde/(1-e^2)^3+...
%     dCDende(2)*cos(f).^2 -...
%     2*C(2)*cos(f).*sin(f).*dfde/(1-e^2)^3+...
%     dCDende(3)*cos(f).^3 -...
%     3*C(3)*cos(f).^2.*sin(f).*dfde/(1-e^2)^3+...
%     dCDende(4)*cos(f).^4 -...
%     4*C(4)*cos(f).^3.*sin(f).*dfde/(1-e^2)^3+...
%     dCDende(5)*cos(f).^5 -...
%     5*C(5)*cos(f).^4.*sin(f).*dfde/(1-e^2)^3)+...
%     (dSDende(1)*sin(f).*cos(f)+...
%     S(1)*(cos(f).^2-sin(f).^2).*dfde/(1-e^2)^3+...
%     dSDende(2)*sin(f).*cos(f).^2+...
%     S(2)*(cos(f).^3-2*sin(f).^2.*cos(f)).*dfde/(1-e^2)^3+...
%     dSDende(3)*sin(f).*cos(f).^3+...
%     S(3)*(cos(f).^4-3*sin(f).^2.*cos(f).^2).*dfde/(1-e^2)^3+...
%     dSDende(4)*sin(f).*cos(f).^4+...
%     S(4)*(cos(f).^5-4*sin(f).^2.*cos(f).^3).*dfde/(1-e^2)^3));

dJ2deT = -3*primary.mu*primary.J2*primary.Re^2/2/a^3/(1-e^2)^3*(1+e*cos(f)).^2.*...
    (dfde.*(-e*(3*sin(aop+f).^2*sin(i)^2-1).*sin(f) + ...
    2*(1+e*cos(f)).*sin(aop+f).*cos(aop+f)*sin(i)^2) +...
    (3*sin(aop+f).^2*sin(i)^2-1).*(2*e*(1+e*cos(f))/(1-e^2) + cos(f)));

% other stuff
[dJ2diF, dJ2diFreq] = dJ2diFour(oe,M,tol,primary);
[dJ2doF, dJ2doFreq] = dJ2doFour(oe,M,tol,primary);
[dJ2dlF, dJ2dlFreq] = dJ2dlFour(oe,M,tol,primary);
[dJ2daF, dJ2daFreq] = dJ2daFour(oe,M,tol,primary);

% Eccentricity
[dJ2deF, dJ2deFreq] = dJ2deFour(oe,M,tol,primary);
 
% Additional verification for da/dt, cuz I don't trust lambda
fRJ2 = -3*primary.mu*primary.J2*primary.Re^2/2/a^4/(1-e^2)^4*(1+e*cos(f)).^4.*...
    (1-3*sin(i)^2*sin(f+aop).^2);
fTHJ2 = -3*primary.mu*primary.J2*primary.Re^2/a^4/(1-e^2)^4*(1+e*cos(f)).^4.*...
    (sin(i)^2*sin(f+aop).*cos(f+aop));
h = sqrt(primary.mu*a*(1-e^2));
dadtForce = 2*a^2/h*e*sin(f).*fRJ2 + 2*a^2/h*(1+e*cos(f)).*fTHJ2;
dadtPot = 2/a/sqrt(primary.mu/a^3)*dJ2dlF;
dadtErr = norm(dadtForce-dadtPot)/length(f);

% Additional verification for dldt
p = a*(1-e^2);
r = p./(1+e*cos(f));
n = sqrt(primary.mu/a^3);
dldtForce = (p*cos(f)-2*r*e)/(n*a^2*e).*fRJ2 -(p+r).*sin(f)/(n*a^2*e).*fTHJ2;
dldtPot = -2/n/a*dJ2daF-(1-e^2)/n/a^2/e*dJ2deF;
dldtErr = norm(dldtForce - dldtPot)/length(f);

% Errors
diErr = norm(dJ2diT-dJ2diF)/length(f)
doErr = norm(dJ2doT-dJ2doF)/length(f)
dlErr = norm(dJ2dlT-dJ2dlF)/length(f)
deErr = norm(dJ2deT-dJ2deF)/length(f)
daErr = norm(dJ2daT-dJ2daF)/length(f)

dJ2daf = 0*M;
dJ2def = 0*M;
dJ2dif = 0*M;
dJ2dof = 0*M;
dJ2dlf = 0*M;
% Check function
for iM = 1:length(M)
    [~, dJ2f] = LpeJ2Fourier([],[oe(1:5),M(iM)],tol);
    dJ2daf(iM) = dJ2f(1);
    dJ2def(iM) = dJ2f(2);
    dJ2dif(iM) = dJ2f(3);
    dJ2dof(iM) = dJ2f(4);
    dJ2dlf(iM) = dJ2f(5);
end

diErrf = norm(dJ2dif-dJ2diF)/length(f)
doErrf = norm(dJ2dof-dJ2doF)/length(f)
dlErrf = norm(dJ2dlf-dJ2dlF)/length(f)
deErrf = norm(dJ2def-dJ2deF)/length(f)
daErrf = norm(dJ2daf-dJ2daF)/length(f)
%% Plot stuff
figure(1)
plot(M,dJ2daf,M,dJ2daF,'--','linewidth',2)
xlabel('M')
xlim([0,2*pi])
xticks([0,pi/2,pi,3*pi/2,2*pi])
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
legend('True','Four')

figure(2)
plot(vecnorm(dJ2daFreq))

figure(3)
plot(dJ2daf-dJ2daF)

%%
function [J2time, J2freq] = J2PotFourier(oe,M,stop,primary)

if nargin < 3 % default tolerance
    tol = 1e-14;
    kMax = Inf;
    nTol = 1e-14;
elseif stop < 1 % stop is tolerance
    tol = stop;
    nTol = stop;
    kMax = Inf;
else            % stop is maximum iterations
    tol = 0;
    kMax = stop;
    nTol = 1e-14;
end
if nargin < 4
    primary = earth();
end
% handle elements vector
a = oe(1);
e = oe(2);
i = oe(3);
raan = oe(4);
aop = oe(5);
b = (1-sqrt(1-e^2))/e;

% constant potential values
R = -primary.mu*primary.J2*primary.Re^2/2/a^3/(1-e^2)^3;
R0 = -1/2*(1-e^2)^(3/2)*(3*cos(i)^2-1);

% constant vectors
m2 = [0:4].';
m3 = [0:6].';
m4 = [0:8].';
m5 = [0:10].';

a2 = 1/2/sqrt(1-e^2)*[1,-4*e,4*e^2+2,-4*e,1].';
a3 = 1/4/(1-e^2)*[1,-6*e,12*e^2+3,-8*e^3-12*e,12*e^2+3,-6*e,1].';
a4 = 1/8/(1-e^2)^(3/2)*[1,-8*e,24*e^2+4,-(32*e^3+24*e),(16*e^4+48*e^2+6),...
    -(32*e^3+24*e),24*e^2+4,-8*e,1].';
a5 = 1/16/(1-e^2)^2*[1,-10*e,40*e^2+5,-(80*e^3+40*e),(80*e^4+120*e^2+10),...
    -(32*e^5+160*e^3+60*e),(80*e^4+120*e^2+10),-(80*e^3+40*e),40*e^2+5,...
    -10*e,1].';

b1 = -1/2*[-1,2*e,0,-2*e,1].';
b2 = -1/4/sqrt(1-e^2)*[-1,4*e,-4*e^2-1,0,4*e^2+1,-4*e,1].';
b3 = -1/8/(1-e^2)*[-1,6*e,-12*e^2-2,8*e^3+6*e,0,-8*e^3-6*e,12*e^2+2,-6*e,1].';
b4 = -1/16/(1-e^2)^(3/2)*[-1,8*e,-24*e^2-3,32*e^3+16*e,-16*e^4-24*e^2-2,0,...
    16*e^4+24*e^2+2,-32*e^3-16*e,24*e^2+3,-8*e,1].';

C1 = (9*e*sin(i)^2*cos(aop)^2-3*e);

C = [9*e^2*sin(i)^2*cos(aop)^2+3*sin(i)^2*(sin(aop)^2-cos(aop)^2)-3*e^2;...
    3*e^3*sin(i)^2*cos(aop)^2+9*e*sin(i)^2*(sin(aop)^2-cos(aop)^2)-e^3;...
    9*e^2*sin(i)^2*(sin(aop)^2-cos(aop)^2);...
    3*e^3*sin(i)^2*(sin(aop)^2-cos(aop)^2)];

S = [6*sin(i)^2*sin(aop)*cos(aop);...
    18*e*sin(i)^2*sin(aop)*cos(aop);...
    18*e^2*sin(i)^2*sin(aop)*cos(aop);...
    6*e^3*sin(i)^2*sin(aop)*cos(aop)];


J2time = zeros(1,length(M));
J2freq = R*[R0;0];

Ak = inf;
Bk = inf;
k = 1;

while max(abs([Ak,Bk])) > tol && k < kMax
    
    n = 0;
    
    g2 = b.^abs(m2+n+k-2);
    g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
    g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
        3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
        3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
    g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
        e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
        5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
        5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);
    
    Ak = C1*2*(1-e^2)/e*besselj(k,k*e) +...
        besselj(n,-k*e)*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
    Bk = besselj(n,-k*e)*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
    dAk = inf;
    dBk = inf;
    n = 1;
    while max(abs([dAk,dBk])) > nTol || n <= k + 5
        g2 = b.^abs(m2+n+k-2);
        g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
        g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
            3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
            3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
        g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
            e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
            5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
            5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);
        
        dAk = besselj(n,-k*e)*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dBk = besselj(n,-k*e)*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        
        g2 = b.^abs(m2-n+k-2);
        g3 = abs(m3-n+k-2).*b.^abs(m3-n+k-3) + e/sqrt(1-e^2)*b.^abs(m3-n+k-2);
        g4 = abs(m4-n+k-3).*abs(m4-n+k-2)/2.*b.^abs(m4-n+k-4) + ...
            3*e*abs(m4-n+k-2)/2/sqrt(1-e^2).*b.^abs(m4-n+k-3) + ...
            3*e^2/2/(1-e^2)*b.^abs(m4-n+k-2);
        g5 = abs(m5-n+k-4).*abs(m5-n+k-3).*abs(m5-n+k-2)/6.*b.^abs(m5-n+k-5) + ...
            e*abs(m5-n+k-3).*abs(m5-n+k-2)/sqrt(1-e^2).*b.^abs(m5-n+k-4) + ...
            5*e^2*abs(m5-n+k-2)/2/(1-e^2).*b.^abs(m5-n+k-3) + ...
            5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5-n+k-2);
        
        dAk = dAk + besselj(-n,-k*e)*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dBk = dBk + besselj(-n,-k*e)*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        
        Ak = Ak + dAk;
        Bk = Bk + dBk;
        n = n+1;
    end
    J2freq = [J2freq R*[Ak;Bk]];
    k = k+1;
end
for k = 0:size(J2freq,2)-1
    J2time = J2time + (J2freq(1,k+1)*cos(k*M) + J2freq(2,k+1)*sin(k*M));
end
end

function [dJ2diTime, dJ2diFreq] = dJ2diFour(oe,M,stop,primary)
% Partial derivative of J2 potential by i - for rate of change of RAAN &
% AOP subject to J2

if nargin < 3 % default tolerance
    tol = 1e-14;
    kMax = Inf;
    nTol = 1e-14;
elseif stop < 1 % stop is tolerance
    tol = stop;
    nTol = stop;
    kMax = Inf;
else            % stop is maximum iterations
    tol = 0;
    kMax = stop;
    nTol = 1e-14;
end
if nargin < 4
    primary = earth();
end
% handle elements vector
a = oe(1);
e = oe(2);
i = oe(3);
raan = oe(4);
aop = oe(5);
b = (1-sqrt(1-e^2))/e;

% constant potential values
R = -primary.mu*primary.J2*primary.Re^2/2/a^3/(1-e^2)^3;
dR0di = 3*(1-e^2)^(3/2)*cos(i)*sin(i);

% constant vectors
m2 = [0:4].';
m3 = [0:6].';
m4 = [0:8].';
m5 = [0:10].';

a2 = 1/2/sqrt(1-e^2)*[1,-4*e,4*e^2+2,-4*e,1].';
a3 = 1/4/(1-e^2)*[1,-6*e,12*e^2+3,-8*e^3-12*e,12*e^2+3,-6*e,1].';
a4 = 1/8/(1-e^2)^(3/2)*[1,-8*e,24*e^2+4,-(32*e^3+24*e),(16*e^4+48*e^2+6),...
    -(32*e^3+24*e),24*e^2+4,-8*e,1].';
a5 = 1/16/(1-e^2)^2*[1,-10*e,40*e^2+5,-(80*e^3+40*e),(80*e^4+120*e^2+10),...
    -(32*e^5+160*e^3+60*e),(80*e^4+120*e^2+10),-(80*e^3+40*e),40*e^2+5,...
    -10*e,1].';

b1 = -1/2*[-1,2*e,0,-2*e,1].';
b2 = -1/4/sqrt(1-e^2)*[-1,4*e,-4*e^2-1,0,4*e^2+1,-4*e,1].';
b3 = -1/8/(1-e^2)*[-1,6*e,-12*e^2-2,8*e^3+6*e,0,-8*e^3-6*e,12*e^2+2,-6*e,1].';
b4 = -1/16/(1-e^2)^(3/2)*[-1,8*e,-24*e^2-3,32*e^3+16*e,-16*e^4-24*e^2-2,0,...
    16*e^4+24*e^2+2,-32*e^3-16*e,24*e^2+3,-8*e,1].';

dC1di = 18*e*sin(i)*cos(i)*cos(aop)^2; % factor multiplying cos(f)

dCdi = [18*e^2*sin(i)*cos(i)*cos(aop)^2+6*sin(i)*cos(i)*(sin(aop)^2-cos(aop)^2);...
    6*e^3*sin(i)*cos(i)*cos(aop)^2+18*e*sin(i)*cos(i)*(sin(aop)^2-cos(aop)^2);...
    18*e^2*sin(i)*cos(i)*(sin(aop)^2-cos(aop)^2);...
    6*e^3*sin(i)*cos(i)*(sin(aop)^2-cos(aop)^2)];

dSdi = [12*sin(i)*cos(i)*sin(aop)*cos(aop);...
    36*e*sin(i)*cos(i)*sin(aop)*cos(aop);...
    36*e^2*sin(i)*cos(i)*sin(aop)*cos(aop);...
    12*e^3*sin(i)*cos(i)*sin(aop)*cos(aop)];


dJ2diTime = zeros(1,length(M));
dJ2diFreq = R*[dR0di;0];

Ak = inf;
Bk = inf;
k = 1;

while max(abs([Ak,Bk])) > tol && k < kMax
    
    n = 0;
    
    g2 = b.^abs(m2+n+k-2);
    g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
    g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
        3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
        3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
    g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
        e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
        5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
        5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);
    
    Ak = dC1di*2*(1-e^2)/e*besselj(k,k*e) +...
        besselj(n,-k*e)*dCdi.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
    Bk = besselj(n,-k*e)*dSdi.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
    dAk = inf;
    dBk = inf;
    n = 1;
    while max(abs([dAk,dBk])) > nTol || n <= k + 5
        g2 = b.^abs(m2+n+k-2);
        g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
        g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
            3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
            3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
        g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
            e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
            5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
            5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);
        
        dAk = besselj(n,-k*e)*dCdi.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dBk = besselj(n,-k*e)*dSdi.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        
        g2 = b.^abs(m2-n+k-2);
        g3 = abs(m3-n+k-2).*b.^abs(m3-n+k-3) + e/sqrt(1-e^2)*b.^abs(m3-n+k-2);
        g4 = abs(m4-n+k-3).*abs(m4-n+k-2)/2.*b.^abs(m4-n+k-4) + ...
            3*e*abs(m4-n+k-2)/2/sqrt(1-e^2).*b.^abs(m4-n+k-3) + ...
            3*e^2/2/(1-e^2)*b.^abs(m4-n+k-2);
        g5 = abs(m5-n+k-4).*abs(m5-n+k-3).*abs(m5-n+k-2)/6.*b.^abs(m5-n+k-5) + ...
            e*abs(m5-n+k-3).*abs(m5-n+k-2)/sqrt(1-e^2).*b.^abs(m5-n+k-4) + ...
            5*e^2*abs(m5-n+k-2)/2/(1-e^2).*b.^abs(m5-n+k-3) + ...
            5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5-n+k-2);
        
        dAk = dAk + besselj(-n,-k*e)*dCdi.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dBk = dBk + besselj(-n,-k*e)*dSdi.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        
        Ak = Ak + dAk;
        Bk = Bk + dBk;
        n = n+1;
    end
    dJ2diFreq = [dJ2diFreq R*[Ak;Bk]];
    k = k+1;
end
for k = 0:size(dJ2diFreq,2)-1
    dJ2diTime = dJ2diTime + (dJ2diFreq(1,k+1)*cos(k*M) + dJ2diFreq(2,k+1)*sin(k*M));
end
end

function [dJ2doTime, dJ2doFreq] = dJ2doFour(oe,M,stop,primary)
% partial derivative of J2 potential by AOP - for rate of change of i & e
% subject to J2

if nargin < 3 % default tolerance
    tol = 1e-14;
    kMax = Inf;
    nTol = 1e-14;
elseif stop < 1 % stop is tolerance
    tol = stop;
    nTol = stop;
    kMax = Inf;
else            % stop is maximum iterations
    tol = 0;
    kMax = stop;
    nTol = 1e-14;
end
if nargin < 4
    primary = earth();
end
% handle elements vector
a = oe(1);
e = oe(2);
i = oe(3);
raan = oe(4);
aop = oe(5);
b = (1-sqrt(1-e^2))/e;

% constant potential values
R = -primary.mu*primary.J2*primary.Re^2/2/a^3/(1-e^2)^3;
dR0do = 0;

% constant vectors
m2 = [0:4].';
m3 = [0:6].';
m4 = [0:8].';
m5 = [0:10].';

a2 = 1/2/sqrt(1-e^2)*[1,-4*e,4*e^2+2,-4*e,1].';
a3 = 1/4/(1-e^2)*[1,-6*e,12*e^2+3,-8*e^3-12*e,12*e^2+3,-6*e,1].';
a4 = 1/8/(1-e^2)^(3/2)*[1,-8*e,24*e^2+4,-(32*e^3+24*e),(16*e^4+48*e^2+6),...
    -(32*e^3+24*e),24*e^2+4,-8*e,1].';
a5 = 1/16/(1-e^2)^2*[1,-10*e,40*e^2+5,-(80*e^3+40*e),(80*e^4+120*e^2+10),...
    -(32*e^5+160*e^3+60*e),(80*e^4+120*e^2+10),-(80*e^3+40*e),40*e^2+5,...
    -10*e,1].';

b1 = -1/2*[-1,2*e,0,-2*e,1].';
b2 = -1/4/sqrt(1-e^2)*[-1,4*e,-4*e^2-1,0,4*e^2+1,-4*e,1].';
b3 = -1/8/(1-e^2)*[-1,6*e,-12*e^2-2,8*e^3+6*e,0,-8*e^3-6*e,12*e^2+2,-6*e,1].';
b4 = -1/16/(1-e^2)^(3/2)*[-1,8*e,-24*e^2-3,32*e^3+16*e,-16*e^4-24*e^2-2,0,...
    16*e^4+24*e^2+2,-32*e^3-16*e,24*e^2+3,-8*e,1].';

dC1do = -18*e*sin(i)^2*sin(aop)*cos(aop); % factor multiplying cos(f)

dCdo = [(-18*e^2+12)*sin(i)^2*sin(aop)*cos(aop);...
    -6*e*(e^2-6)*sin(i)^2*sin(aop)*cos(aop);...
    36*e^2*sin(i)^2*sin(aop)*cos(aop);...
    12*e^3*sin(i)^2*sin(aop)*cos(aop)];

dSdo = [12*sin(i)^2*cos(aop)^2-6*sin(i)^2;...
    18*e*sin(i)^2*(2*cos(aop)^2-1);...
    18*e^2*sin(i)^2*(2*cos(aop)^2-1);...
    6*e^3*sin(i)^2*(2*cos(aop)^2-1)];


dJ2doTime = zeros(1,length(M));
dJ2doFreq = R*[dR0do;0];

Ak = inf;
Bk = inf;
k = 1;

while max(abs([Ak,Bk])) > tol && k < kMax
    
    n = 0;
    
    g2 = b.^abs(m2+n+k-2);
    g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
    g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
        3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
        3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
    g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
        e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
        5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
        5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);
    
    Ak = dC1do*2*(1-e^2)/e*besselj(k,k*e) +...
        besselj(n,-k*e)*dCdo.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
    Bk = besselj(n,-k*e)*dSdo.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
    dAk = inf;
    dBk = inf;
    n = 1;
    while max(abs([dAk,dBk])) > nTol || n <= k + 5
        g2 = b.^abs(m2+n+k-2);
        g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
        g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
            3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
            3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
        g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
            e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
            5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
            5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);
        
        dAk = besselj(n,-k*e)*dCdo.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dBk = besselj(n,-k*e)*dSdo.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        
        g2 = b.^abs(m2-n+k-2);
        g3 = abs(m3-n+k-2).*b.^abs(m3-n+k-3) + e/sqrt(1-e^2)*b.^abs(m3-n+k-2);
        g4 = abs(m4-n+k-3).*abs(m4-n+k-2)/2.*b.^abs(m4-n+k-4) + ...
            3*e*abs(m4-n+k-2)/2/sqrt(1-e^2).*b.^abs(m4-n+k-3) + ...
            3*e^2/2/(1-e^2)*b.^abs(m4-n+k-2);
        g5 = abs(m5-n+k-4).*abs(m5-n+k-3).*abs(m5-n+k-2)/6.*b.^abs(m5-n+k-5) + ...
            e*abs(m5-n+k-3).*abs(m5-n+k-2)/sqrt(1-e^2).*b.^abs(m5-n+k-4) + ...
            5*e^2*abs(m5-n+k-2)/2/(1-e^2).*b.^abs(m5-n+k-3) + ...
            5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5-n+k-2);
        
        dAk = dAk + besselj(-n,-k*e)*dCdo.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dBk = dBk + besselj(-n,-k*e)*dSdo.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        
        Ak = Ak + dAk;
        Bk = Bk + dBk;
        n = n+1;
    end
    dJ2doFreq = [dJ2doFreq R*[Ak;Bk]];
    k = k+1;
end
for k = 0:size(dJ2doFreq,2)-1
    dJ2doTime = dJ2doTime + (dJ2doFreq(1,k+1)*cos(k*M) + dJ2doFreq(2,k+1)*sin(k*M));
end
end

function [dJ2dlTime, dJ2dlFreq] = dJ2dlFour(oe,M,stop,primary)

if nargin < 3 % default tolerance
    tol = 1e-14;
    kMax = Inf;
    nTol = 1e-14;
elseif stop < 1 % stop is tolerance
    tol = stop;
    nTol = stop;
    kMax = Inf;
else            % stop is maximum iterations
    tol = 0;
    kMax = stop;
    nTol = 1e-14;
end
if nargin < 4
    primary = earth();
end
% handle elements vector
a = oe(1);
e = oe(2);
i = oe(3);
raan = oe(4);
aop = oe(5);
b = (1-sqrt(1-e^2))/e;

% constant potential values
R = -primary.mu*primary.J2*primary.Re^2/2/a^3/(1-e^2)^3;
dR0dl = 0;

% constant vectors
m2 = [0:4].';
m3 = [0:6].';
m4 = [0:8].';
m5 = [0:10].';

a2 = 1/2/sqrt(1-e^2)*[1,-4*e,4*e^2+2,-4*e,1].';
a3 = 1/4/(1-e^2)*[1,-6*e,12*e^2+3,-8*e^3-12*e,12*e^2+3,-6*e,1].';
a4 = 1/8/(1-e^2)^(3/2)*[1,-8*e,24*e^2+4,-(32*e^3+24*e),(16*e^4+48*e^2+6),...
    -(32*e^3+24*e),24*e^2+4,-8*e,1].';
a5 = 1/16/(1-e^2)^2*[1,-10*e,40*e^2+5,-(80*e^3+40*e),(80*e^4+120*e^2+10),...
    -(32*e^5+160*e^3+60*e),(80*e^4+120*e^2+10),-(80*e^3+40*e),40*e^2+5,...
    -10*e,1].';

b1 = -1/2*[-1,2*e,0,-2*e,1].';
b2 = -1/4/sqrt(1-e^2)*[-1,4*e,-4*e^2-1,0,4*e^2+1,-4*e,1].';
b3 = -1/8/(1-e^2)*[-1,6*e,-12*e^2-2,8*e^3+6*e,0,-8*e^3-6*e,12*e^2+2,-6*e,1].';
b4 = -1/16/(1-e^2)^(3/2)*[-1,8*e,-24*e^2-3,32*e^3+16*e,-16*e^4-24*e^2-2,0,...
    16*e^4+24*e^2+2,-32*e^3-16*e,24*e^2+3,-8*e,1].';

C1 = (9*e*sin(i)^2*cos(aop)^2-3*e);

C = [9*e^2*sin(i)^2*cos(aop)^2+3*sin(i)^2*(sin(aop)^2-cos(aop)^2)-3*e^2;...
    3*e^3*sin(i)^2*cos(aop)^2+9*e*sin(i)^2*(sin(aop)^2-cos(aop)^2)-e^3;...
    9*e^2*sin(i)^2*(sin(aop)^2-cos(aop)^2);...
    3*e^3*sin(i)^2*(sin(aop)^2-cos(aop)^2)];

S = [6*sin(i)^2*sin(aop)*cos(aop);...
    18*e*sin(i)^2*sin(aop)*cos(aop);...
    18*e^2*sin(i)^2*sin(aop)*cos(aop);...
    6*e^3*sin(i)^2*sin(aop)*cos(aop)];


dJ2dlTime = zeros(1,length(M));
dJ2dlFreq = R*[dR0dl;0];

Ak = inf;
Bk = inf;
k = 1;

while max(abs([Ak,Bk])) > tol && k < kMax
    
    n = 0;
    
    g2 = b.^abs(m2+n+k-2);
    g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
    g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
        3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
        3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
    g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
        e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
        5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
        5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);
    
    Ak = C1*2*(1-e^2)/e*besselj(k,k*e) +...
        besselj(n,-k*e)*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
    Bk = besselj(n,-k*e)*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
    dAk = inf;
    dBk = inf;
    n = 1;
    while max(abs([dAk,dBk])) > nTol || n <= k + 5
        g2 = b.^abs(m2+n+k-2);
        g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
        g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
            3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
            3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
        g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
            e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
            5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
            5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);
        
        dAk = besselj(n,-k*e)*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dBk = besselj(n,-k*e)*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        
        g2 = b.^abs(m2-n+k-2);
        g3 = abs(m3-n+k-2).*b.^abs(m3-n+k-3) + e/sqrt(1-e^2)*b.^abs(m3-n+k-2);
        g4 = abs(m4-n+k-3).*abs(m4-n+k-2)/2.*b.^abs(m4-n+k-4) + ...
            3*e*abs(m4-n+k-2)/2/sqrt(1-e^2).*b.^abs(m4-n+k-3) + ...
            3*e^2/2/(1-e^2)*b.^abs(m4-n+k-2);
        g5 = abs(m5-n+k-4).*abs(m5-n+k-3).*abs(m5-n+k-2)/6.*b.^abs(m5-n+k-5) + ...
            e*abs(m5-n+k-3).*abs(m5-n+k-2)/sqrt(1-e^2).*b.^abs(m5-n+k-4) + ...
            5*e^2*abs(m5-n+k-2)/2/(1-e^2).*b.^abs(m5-n+k-3) + ...
            5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5-n+k-2);
        
        dAk = dAk + besselj(-n,-k*e)*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dBk = dBk + besselj(-n,-k*e)*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        
        Ak = Ak + dAk;
        Bk = Bk + dBk;
        n = n+1;
    end
    dJ2dlFreq = [dJ2dlFreq R*[Ak;Bk]];
    k = k+1;
end
for k = 0:size(dJ2dlFreq,2)-1
    dJ2dlTime = dJ2dlTime + (-k*dJ2dlFreq(1,k+1)*sin(k*M) + k*dJ2dlFreq(2,k+1)*cos(k*M));
end
end

function [dJ2detime, dJ2deFreq] = dJ2deFour(oe,M,stop,primary)

if nargin < 3 % default tolerance
    tol = 1e-14;
    kMax = Inf;
    nTol = 1e-14;
elseif stop < 1 % stop is tolerance
    tol = stop;
    nTol = stop;
    kMax = Inf;
else            % stop is maximum iterations
    tol = 0;
    kMax = stop;
    nTol = 1e-14;
end
if nargin < 4
    primary = earth();
end
% handle elements vector
a = oe(1);
e = oe(2);
i = oe(3);
raan = oe(4);
aop = oe(5);
b = (1-sqrt(1-e^2))/e;
dbde = (1+b^2)/2/sqrt(1-e^2);

% constant potential values
R0Tilde = -primary.mu*primary.J2*primary.Re^2/2/a^3;
dRBarde = -3*e*(3*cos(i)^2-1)/2/(1-e^2)^(5/2); %

% constant vectors
m2 = [0:4].';
m3 = [0:6].';
m4 = [0:8].';
m5 = [0:10].';

a2 = [1,-4*e,4*e^2+2,-4*e,1].';
a3 = [1,-6*e,12*e^2+3,-8*e^3-12*e,12*e^2+3,-6*e,1].';
a4 = [1,-8*e,24*e^2+4,-(32*e^3+24*e),(16*e^4+48*e^2+6),...
    -(32*e^3+24*e),24*e^2+4,-8*e,1].';
a5 = [1,-10*e,40*e^2+5,-(80*e^3+40*e),(80*e^4+120*e^2+10),...
    -(32*e^5+160*e^3+60*e),(80*e^4+120*e^2+10),-(80*e^3+40*e),40*e^2+5,...
    -10*e,1].';

b1 = [-1,2*e,0,-2*e,1].';
b2 = [-1,4*e,-4*e^2-1,0,4*e^2+1,-4*e,1].';
b3 = [-1,6*e,-12*e^2-2,8*e^3+6*e,0,-8*e^3-6*e,12*e^2+2,-6*e,1].';
b4 = [-1,8*e,-24*e^2-3,32*e^3+16*e,-16*e^4-24*e^2-2,0,16*e^4+24*e^2+2,...
    -32*e^3-16*e,24*e^2+3,-8*e,1].';

da2de = [0,-4,8*e,-4,0].';
da3de = [0,-6,24*e,(-24*e^2-12),24*e,-6,0].';
da4de = [0,-8,48*e,(-96*e^2-24),(64*e^3+96*e),(-96*e^2-24),48*e,-8,0].';
da5de = [0,-10,80*e,(-240*e^2-40),(320*e^3+240*e),(-160*e^4-480*e^2-60),...
        (320*e^3+240*e),(-240*e^2-40),80*e,-10,0].';
    
db1de = [0,2,0,-2,0].';
db2de = [0,4,-8*e,0,8*e,-4,0].';
db3de = [0,6,-24*e,24*e^2+6,0,-24*e^2-6,24*e,-6,0].';
db4de = [0,8,-48*e,96*e^2+16,-64*e^3-48*e,0,64*e^3+48*e,-96*e^2-16,48*e,-8,0].';

    
C1xA1 = 6*(3*sin(i)^2*cos(aop)^2-1)/(1-e^2)^2;

dC1xA1de =  24*e*(3*sin(i)^2*cos(aop)^2-1)/(1-e^2)^3;

% C & S incorporate denominators
CC = [(9*e^2*sin(i)^2*cos(aop)^2+3*sin(i)^2*(sin(aop)^2-cos(aop)^2)-3*e^2)/2/(1-e^2)^(7/2);
    (3*e^3*sin(i)^2*cos(aop)^2+9*e*sin(i)^2*(sin(aop)^2-cos(aop)^2)-e^3)/4/(1-e^2)^4;
    (9*e^2*sin(i)^2*(sin(aop)^2-cos(aop)^2))/8/(1-e^2)^(9/2);
    (3*e^3*sin(i)^2*(sin(aop)^2-cos(aop)^2))/16/(1-e^2)^5];

SS = -[6*sin(i)^2*sin(aop)*cos(aop)/2/(1-e^2)^3;
    18*e*sin(i)^2*sin(aop)*cos(aop)/4/(1-e^2)^(7/2);
    18*e^2*sin(i)^2*sin(aop)*cos(aop)/8/(1-e^2)^4;
    6*e^3*sin(i)^2*sin(aop)*cos(aop)/16/(1-e^2)^(9/2)];

% Derivative of C & S including denominators
dCCde = [3*e*(15*e^2*sin(i)^2*cos(aop)^2 +7*sin(i)^2*sin(aop)^2 -sin(i)^2*cos(aop)^2 -5*e^2 -2)/2/(1-e^2)^(9/2);
        ((15*e^4-54*e^2-9)*sin(i)^2*cos(aop)^2 + (63*e^2+9)*sin(i)^2*sin(aop)^2 -5*e^4 -3*e^2)/4/(1-e^2)^5;
        -9*e*sin(i)^2*(cos(aop)^2-sin(aop)^2)*(7*e^2+2)/8/(1-e^2)^(11/2);
        -3*e^2*sin(i)^2*(cos(aop)^2-sin(aop)^2)*(7*e^2+3)/16/(1-e^2)^6];

dSSde = [-18*e*sin(i)^2*sin(aop)*cos(aop)/(1-e^2)^4;
         -9*sin(i)^2*sin(aop)*cos(aop)*(6*e^2+1)/2/(1-e^2)^(9/2);
         -9*e*sin(i)^2*sin(aop)*cos(aop)*(3*e^2+1)/2/(1-e^2)^5;
         -9*e^2*sin(i)^2*sin(aop)*cos(aop)*(2*e^2+1)/8/(1-e^2)^(11/2)];


dJ2detime = zeros(1,length(M));
dJ2deFreq = R0Tilde*[dRBarde;0];

Ak = inf;
Bk = inf;
k = 1;

while max(abs([Ak,Bk])) > tol && k < kMax
    
    n = 0;
    
    g2 = b.^abs(m2+n+k-2);
    g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
    g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
        3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
        3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
    g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
        e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
        5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
        5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);
    
    dg2de = abs(m2+n+k-2).*b.^(abs(m2+n+k-2))/e/sqrt(1-e^2);
    dg3de = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^(abs(m3+n+k-3)) + ...
            e/sqrt(1-e^2)*b.^(abs(m3+n+k-2)))/e/sqrt(1-e^2) + ...
            b.^abs(m3+n+k-2)/(1-e^2)^(3/2);
    dg4de = abs(m4+n+k-2).*(3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2) +...
            abs(m4+n+k-3).*(3*e/2/sqrt(1-e^2)*b.^abs(m4+n+k-3) +...
            abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))/e/sqrt(1-e^2) +...
            3/2*abs(m4+n+k-2)/(1-e^2)^(3/2).*b.^abs(m4+n+k-3) +...
            3*e/(1-e^2)^2*b.^abs(m4+n+k-2);
    dg5de = abs(m5+n+k-2).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(m5+n+k-2) + ...
            abs(m5+n+k-3).*(5*e^2/2/(1-e^2)*b.^abs(m5+n+k-3) + ...
            abs(m5+n+k-4).*(e/sqrt(1-e^2)*b.^abs(m5+n+k-4) + ...
            abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))/e/sqrt(1-e^2) + ...
            abs(m5+n+k-3).*abs(m5+n+k-2)/(1-e^2)^(3/2).*b.^abs(m5+n+k-4) +...
            5*abs(m5+n+k-2)*e/(1-e^2)^2.*b.^abs(m5+n+k-3) +...
            15*e^2/2/(1-e^2)^(5/2)*b.^abs(m5+n+k-2);
    
    
    dJkde = 0.5*(besselj(k-1,k*e) - besselj(k+1,k*e));
    Jk = besselj(k,k*e);
    
    Jn = besselj(n,-k*e);
    dJnde = 0.5*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
    
    Ak = (C1xA1*k*dJkde + dC1xA1de*Jk) + ...
        -k*dJnde*CC.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]+...
        Jn*(dCCde.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] + ...
        CC.'*[da2de.'*g2+a2.'*dg2de; da3de.'*g3+a3.'*dg3de; da4de.'*g4+a4.'*dg4de; da5de.'*g5+a5.'*dg5de]);
    
    Bk = -k*dJnde*SS.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5]+...
        Jn*(dSSde.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] + ...
        SS.'*[db1de.'*g2+b1.'*dg2de; db2de.'*g3+b2.'*dg3de; db3de.'*g4+b3.'*dg4de; db4de.'*g5+b4.'*dg5de]);
    
    dAk = inf;
    dBk = inf;
    n = 1;
    while max(abs([dAk,dBk])) > nTol || n <= k + 5
        % Positive n
        g2 = b.^abs(m2+n+k-2);
        g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
        g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
            3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
            3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
        g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
            e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
            5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
            5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);
        
        dg2de = abs(m2+n+k-2).*b.^abs(m2+n+k-2)/e/sqrt(1-e^2);
        dg3de = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^abs(m3+n+k-3) + ...
                e/sqrt(1-e^2)*b.^abs(m3+n+k-2))/e/sqrt(1-e^2) + ...
                b.^abs(m3+n+k-2)/(1-e^2)^(3/2);
        dg4de = abs(m4+n+k-2).*(3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2) +...
                abs(m4+n+k-3).*(3*e/2/sqrt(1-e^2)*b.^abs(m4+n+k-3) +...
                abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))/e/sqrt(1-e^2) +...
                3/2*abs(m4+n+k-2)/(1-e^2)^(3/2).*b.^abs(m4+n+k-3) +...
                3*e/(1-e^2)^2*b.^abs(m4+n+k-2);
        dg5de = abs(m5+n+k-2).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(m5+n+k-2) + ...
                abs(m5+n+k-3).*(5*e^2/2/(1-e^2)*b.^abs(m5+n+k-3) + ...
                abs(m5+n+k-4).*(e/sqrt(1-e^2)*b.^abs(m5+n+k-4) + ...
                abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))/e/sqrt(1-e^2) + ...
                abs(m5+n+k-3).*abs(m5+n+k-2)/(1-e^2)^(3/2).*b.^abs(m5+n+k-4) +...
                5*abs(m5+n+k-2)*e/(1-e^2)^2.*b.^abs(m5+n+k-3) +...
                15*e^2/2/(1-e^2)^(5/2)*b.^abs(m5+n+k-2);
        
        Jn = besselj(n,-k*e);
        dJnde = 0.5*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
        
        
        dAk = -k*dJnde*CC.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]+...
        Jn*(dCCde.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] + ...
        CC.'*[da2de.'*g2+a2.'*dg2de; da3de.'*g3+a3.'*dg3de; da4de.'*g4+a4.'*dg4de; da5de.'*g5+a5.'*dg5de]);
        
        dBk = -k*dJnde*SS.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5]+...
        Jn*(dSSde.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] + ...
        SS.'*[db1de.'*g2+b1.'*dg2de; db2de.'*g3+b2.'*dg3de; db3de.'*g4+b3.'*dg4de; db4de.'*g5+b4.'*dg5de]);
        
        % negative n
        g2 = b.^abs(m2-n+k-2);
        g3 = abs(m3-n+k-2).*b.^abs(m3-n+k-3) + e/sqrt(1-e^2)*b.^abs(m3-n+k-2);
        g4 = abs(m4-n+k-3).*abs(m4-n+k-2)/2.*b.^abs(m4-n+k-4) + ...
            3*e*abs(m4-n+k-2)/2/sqrt(1-e^2).*b.^abs(m4-n+k-3) + ...
            3*e^2/2/(1-e^2)*b.^abs(m4-n+k-2);
        g5 = abs(m5-n+k-4).*abs(m5-n+k-3).*abs(m5-n+k-2)/6.*b.^abs(m5-n+k-5) + ...
            e*abs(m5-n+k-3).*abs(m5-n+k-2)/sqrt(1-e^2).*b.^abs(m5-n+k-4) + ...
            5*e^2*abs(m5-n+k-2)/2/(1-e^2).*b.^abs(m5-n+k-3) + ...
            5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5-n+k-2);
        
        dg2de = abs(m2-n+k-2).*b.^(abs(m2-n+k-2))/e/sqrt(1-e^2);
        dg3de = abs(m3-n+k-2).*(abs(m3-n+k-3).*b.^(abs(m3-n+k-3)) + ...
                e/sqrt(1-e^2)*b.^(abs(m3-n+k-2)))/e/sqrt(1-e^2) + ...
                b.^abs(m3-n+k-2)/(1-e^2)^(3/2);
        dg4de = abs(m4-n+k-2).*(3*e^2/2/(1-e^2)*b.^(abs(m4-n+k-2)) +...
                abs(m4-n+k-3).*(3*e/2/sqrt(1-e^2)*b.^(abs(m4-n+k-3)) +...
                abs(m4-n+k-4)/2.*b.^abs(m4-n+k-4)))/e/sqrt(1-e^2) +...
                3/2*abs(m4-n+k-2)/(1-e^2)^(3/2).*b.^abs(m4-n+k-3) +...
                3*e/(1-e^2)^2*b.^abs(m4-n+k-2);
        dg5de = abs(m5-n+k-2).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(m5-n+k-2) + ...
                abs(m5-n+k-3).*(5*e^2/2/(1-e^2)*b.^abs(m5-n+k-3) + ...
                abs(m5-n+k-4).*(e/sqrt(1-e^2)*b.^abs(m5-n+k-4) + ...
                abs(m5-n+k-5)/6.*b.^abs(m5-n+k-5))))/e/sqrt(1-e^2) + ...
                abs(m5-n+k-3).*abs(m5-n+k-2)/(1-e^2)^(3/2).*b.^abs(m5-n+k-4) +...
                5*abs(m5-n+k-2)*e/(1-e^2)^2.*b.^abs(m5-n+k-3) +...
                15*e^2/2/(1-e^2)^(5/2)*b.^abs(m5-n+k-2);
            
        Jn = besselj(n,k*e);
        dJnde = 0.5*(besselj(n-1,k*e) - besselj(n+1,k*e));
        
        dAk = dAk + k*dJnde*CC.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]+...
        Jn*(dCCde.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] + ...
        CC.'*[da2de.'*g2+a2.'*dg2de; da3de.'*g3+a3.'*dg3de; da4de.'*g4+a4.'*dg4de; da5de.'*g5+a5.'*dg5de]);
    
        dBk = dBk + k*dJnde*SS.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5]+...
        Jn*(dSSde.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] + ...
        SS.'*[db1de.'*g2+b1.'*dg2de; db2de.'*g3+b2.'*dg3de; db3de.'*g4+b3.'*dg4de; db4de.'*g5+b4.'*dg5de]);
        
        Ak = Ak + dAk;
        Bk = Bk + dBk;
        n = n+1;
    end
    dJ2deFreq = [dJ2deFreq R0Tilde*[Ak;Bk]]; %#ok<AGROW> 
    k = k+1;
end
for k = 0:size(dJ2deFreq,2)-1
    dJ2detime = dJ2detime + (dJ2deFreq(1,k+1)*cos(k*M) + dJ2deFreq(2,k+1)*sin(k*M));
end
end

function [dJ2datime, dJ2dafreq] = dJ2daFour(oe,M,stop,primary)

if nargin < 3 % default tolerance
    tol = 1e-14;
    kMax = Inf;
    nTol = 1e-14;
elseif stop < 1 % stop is tolerance
    tol = stop;
    nTol = stop;
    kMax = Inf;
else            % stop is maximum iterations
    tol = 0;
    kMax = stop;
    nTol = 1e-14;
end
if nargin < 4
    primary = earth();
end
% handle elements vector
a = oe(1);
e = oe(2);
i = oe(3);
raan = oe(4);
aop = oe(5);
b = (1-sqrt(1-e^2))/e;

% constant potential values
R = 3*primary.mu*primary.J2*primary.Re^2/2/a^4/(1-e^2)^3;
R0 = -1/2*(1-e^2)^(3/2)*(3*cos(i)^2-1);

% constant vectors
m2 = [0:4].';
m3 = [0:6].';
m4 = [0:8].';
m5 = [0:10].';

a2 = 1/2/sqrt(1-e^2)*[1,-4*e,4*e^2+2,-4*e,1].';
a3 = 1/4/(1-e^2)*[1,-6*e,12*e^2+3,-8*e^3-12*e,12*e^2+3,-6*e,1].';
a4 = 1/8/(1-e^2)^(3/2)*[1,-8*e,24*e^2+4,-(32*e^3+24*e),(16*e^4+48*e^2+6),...
    -(32*e^3+24*e),24*e^2+4,-8*e,1].';
a5 = 1/16/(1-e^2)^2*[1,-10*e,40*e^2+5,-(80*e^3+40*e),(80*e^4+120*e^2+10),...
    -(32*e^5+160*e^3+60*e),(80*e^4+120*e^2+10),-(80*e^3+40*e),40*e^2+5,...
    -10*e,1].';

b1 = -1/2*[-1,2*e,0,-2*e,1].';
b2 = -1/4/sqrt(1-e^2)*[-1,4*e,-4*e^2-1,0,4*e^2+1,-4*e,1].';
b3 = -1/8/(1-e^2)*[-1,6*e,-12*e^2-2,8*e^3+6*e,0,-8*e^3-6*e,12*e^2+2,-6*e,1].';
b4 = -1/16/(1-e^2)^(3/2)*[-1,8*e,-24*e^2-3,32*e^3+16*e,-16*e^4-24*e^2-2,0,...
    16*e^4+24*e^2+2,-32*e^3-16*e,24*e^2+3,-8*e,1].';

C1 = (9*e*sin(i)^2*cos(aop)^2-3*e);

C = [9*e^2*sin(i)^2*cos(aop)^2+3*sin(i)^2*(sin(aop)^2-cos(aop)^2)-3*e^2;...
    3*e^3*sin(i)^2*cos(aop)^2+9*e*sin(i)^2*(sin(aop)^2-cos(aop)^2)-e^3;...
    9*e^2*sin(i)^2*(sin(aop)^2-cos(aop)^2);...
    3*e^3*sin(i)^2*(sin(aop)^2-cos(aop)^2)];

S = [6*sin(i)^2*sin(aop)*cos(aop);...
    18*e*sin(i)^2*sin(aop)*cos(aop);...
    18*e^2*sin(i)^2*sin(aop)*cos(aop);...
    6*e^3*sin(i)^2*sin(aop)*cos(aop)];


dJ2datime = zeros(1,length(M));
dJ2dafreq = R*[R0;0];

Ak = inf;
Bk = inf;
k = 1;

while max(abs([Ak,Bk])) > tol && k < kMax
    
    n = 0;
    
    g2 = b.^abs(m2+n+k-2);
    g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
    g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
        3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
        3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
    g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
        e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
        5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
        5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);
    
    Ak = C1*2*(1-e^2)/e*besselj(k,k*e) +...
        besselj(n,-k*e)*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
    Bk = besselj(n,-k*e)*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
    dAk = inf;
    dBk = inf;
    n = 1;
    while max(abs([dAk,dBk])) > nTol || n <= k + 5
        g2 = b.^abs(m2+n+k-2);
        g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
        g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
            3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
            3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
        g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
            e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
            5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
            5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);
        
        dAk = besselj(n,-k*e)*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dBk = besselj(n,-k*e)*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        
        g2 = b.^abs(m2-n+k-2);
        g3 = abs(m3-n+k-2).*b.^abs(m3-n+k-3) + e/sqrt(1-e^2)*b.^abs(m3-n+k-2);
        g4 = abs(m4-n+k-3).*abs(m4-n+k-2)/2.*b.^abs(m4-n+k-4) + ...
            3*e*abs(m4-n+k-2)/2/sqrt(1-e^2).*b.^abs(m4-n+k-3) + ...
            3*e^2/2/(1-e^2)*b.^abs(m4-n+k-2);
        g5 = abs(m5-n+k-4).*abs(m5-n+k-3).*abs(m5-n+k-2)/6.*b.^abs(m5-n+k-5) + ...
            e*abs(m5-n+k-3).*abs(m5-n+k-2)/sqrt(1-e^2).*b.^abs(m5-n+k-4) + ...
            5*e^2*abs(m5-n+k-2)/2/(1-e^2).*b.^abs(m5-n+k-3) + ...
            5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5-n+k-2);
        
        dAk = dAk + besselj(-n,-k*e)*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dBk = dBk + besselj(-n,-k*e)*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        
        Ak = Ak + dAk;
        Bk = Bk + dBk;
        n = n+1;
    end
    dJ2dafreq = [dJ2dafreq R*[Ak;Bk]];
    k = k+1;
end
for k = 0:size(dJ2dafreq,2)-1
    dJ2datime = dJ2datime + (dJ2dafreq(1,k+1)*cos(k*M) + dJ2dafreq(2,k+1)*sin(k*M));
end
end