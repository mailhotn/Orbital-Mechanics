tol = 1e-14;
primary = earth();
%% Orbital Elements
oe = [8000, 0.3, 30*pi/180, 0, 80*pi/180, 0];
a = oe(1);
e = oe(2);
i = oe(3);
raan = oe(4);
aop = oe(5);
M = 0:0.001:2*pi;
f = pi/180*me2ta(180/pi*M,e,tol);
%% Calculate Potentials

[J2F,J2freq] = J2PotFourier(oe,M,tol);

J2T = -primary.mu*primary.J2*primary.Re^2/2/a^3/(1-e^2)^3*(1+e*cos(f)).^3.*...
    (3*sin(aop+f).^2*sin(i)^2-1);

err = norm(J2T-J2F)/length(f)
%% Plot stuff
figure(1)
plot(M,J2T,M,J2F,'--','linewidth',2)

figure(2)
plot(vecnorm(J2freq))

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
% R*R0 = mean of R - constant factor
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
    % Initial value is of n=0, and C1*A1
    Ak = (9*e*sin(i)^2*cos(aop)^2-3*e)*2*(1-e^2)/e*besselj(k,k*e) +...
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
