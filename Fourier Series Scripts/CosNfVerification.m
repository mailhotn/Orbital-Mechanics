clear
e = 0.1;
M = 0:0.001:2*pi;
tol = 1e-14;
%% Calculate N-R, Battin

fF = fFourier(M,e,tol);

tic
fT = pi/180*me2ta(180/pi*M,e,tol);
c2T = cos(fT).^2;
c3T = cos(fT).^3;
c4T = cos(fT).^4;
c5T = cos(fT).^5;
c1s1T = cos(fT).*sin(fT);
c2s1T = cos(fT).^2.*sin(fT);
c3s1T = cos(fT).^3.*sin(fT);
c4s1T = cos(fT).^4.*sin(fT);
toc

tic
[c2F, BkVec2] = Cos2f(M,e,tol);
toc
tic
[c3F, BkVec3] = Cos3f(M,e,tol);
toc
tic
[c4F, BkVec4] = Cos4f(M,e,tol);
toc
tic
[c5F, BkVec5] = Cos5f(M,e,tol);
toc

tic
[c1s1F,BkVec11] = CosfSinf(M,e,tol);
toc
tic
[c2s1F,BkVec21] = Cos2fSinf(M,e,tol);
toc
tic
[c3s1F,BkVec31] = Cos3fSinf(M,e,tol);
toc
tic
[c4s1F,BkVec41] = Cos4fSinf(M,e,tol);
toc

errorc2f = norm(c2T-c2F)/length(c2T)
errorc3f = norm(c3T-c3F)/length(c3T)
errorc4f = norm(c4T-c4F)/length(c4T)
errorc5f = norm(c5T-c5F)/length(c5T)
errorc1s1F = norm(c1s1T-c1s1F)/length(c1s1T)
errorc2s1F = norm(c2s1T-c2s1F)/length(c2s1T)
errorc3s1F = norm(c3s1T-c3s1F)/length(c3s1T)
errorc4s1F = norm(c4s1T-c4s1F)/length(c4s1T)
%% Plot Stuff
figure(1) % M & E
plot(M,fT,M,fF,'--','linewidth',2)
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
yticks(0:pi/2:2*pi)
yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
axis([0 2*pi 0 2*pi])

figure(2)
plot(M,c2T,M,c2F,'--','linewidth',2)
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('M')
axis([0 2*pi -1 1])
legend('N-R','Fourier','location','best')

figure(3)
plot(M,c3T,M,c3F,'--','linewidth',2)
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('M')
axis([0 2*pi -1 1])
legend('N-R','Fourier','location','best')

figure(4)
plot(M,c4T,M,c4F,'--','linewidth',2)
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('M')
axis([0 2*pi -1 1])
legend('N-R','Fourier','location','best')

figure(5)
plot(M,c5T,M,c5F,'--','linewidth',2)
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('M')
axis([0 2*pi -1 1])
legend('N-R','Fourier','location','best')

figure(6)
plot(M,c1s1T,M,c1s1F,'--','linewidth',2)
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('M')
axis([0 2*pi -1 1])
legend('N-R','Fourier','location','best')

figure(7)
plot(M,c2s1T,M,c2s1F,'--','linewidth',2)
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('M')
axis([0 2*pi -1 1])
legend('N-R','Fourier','location','best')

figure(8)
plot(M,c3s1T,M,c3s1F,'--','linewidth',2)
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('M')
axis([0 2*pi -1 1])
legend('N-R','Fourier','location','best')

figure(9)
plot(M,c4s1T,M,c4s1F,'--','linewidth',2)
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('M')
axis([0 2*pi -1 1])
legend('N-R','Fourier','location','best')

figure(10)
plot(0:length(BkVec2)-1,BkVec2,'.',...
    0:length(BkVec3)-1,BkVec3,'.',...
    0:length(BkVec4)-1,BkVec4,'.',...
    0:length(BkVec5)-1,BkVec5,'.',...
    1:length(BkVec11),BkVec11,'o',...
    1:length(BkVec21),BkVec21,'o',...
    1:length(BkVec31),BkVec31,'o',...
    1:length(BkVec41),BkVec41,'o')
legend('c2','c3','c4','c5','c1s1','c2s1','c3s1','c4s1')
%%
function [c2f,BkVec] = Cos2f(M,e,stop)
if nargin < 3 % default tolerance
    tol = 1e-14;
    tol = nTol;
    kMax = Inf;
elseif stop < 1 % stop is tolerance
    tol = stop;
    nTol = stop;
    kMax = Inf;
else            % stop is maximum iterations
    tol = 0;
    kMax = stop;
    nTol = 1e-14;
end
b = (1-sqrt(1-e^2))/e;
c2f = zeros(1,length(M));
B0 = b*(2*e^2-4*e*b+b^2+1)/e/(1-b^2); % constant term
BkVec = [B0];

Bk = inf;
k = 1;
while abs(Bk) > tol && k < kMax% Iterate until Bk is small enough

    dBk = inf;
    n = 1;
    Bk = 1/2/sqrt(1-e^2)*besselj(0,-k*e)*(b^abs(k+2)-4*e*b^abs(k+1)+...
        (4*e^2+2)*b^abs(k)-4*e*b^abs(k-1)+b^abs(k-2));
    while abs(dBk) > nTol || n <= k + 2% Iterate until Bk converges
        dBk = 1/2/sqrt(1-e^2)*(besselj(n,-k*e)*(b^abs(n+k+2)-4*e*b^abs(n+k+1)+...
            (4*e^2+2)*b^abs(n+k)-4*e*b^abs(n+k-1)+b^abs(n+k-2)) + ...
            besselj(-n,-k*e)*(b^abs(-n+k+2)-4*e*b^abs(-n+k+1)+...
            (4*e^2+2)*b^abs(-n+k)-4*e*b^abs(-n+k-1)+b^abs(-n+k-2)));
        
        Bk = Bk + dBk;
        n = n+1;
    end
    BkVec = [BkVec; Bk];
    k = k+1;
end
for k = 0:length(BkVec)-1
    c2f = c2f + BkVec(k+1)*cos(k*M);
end
end

function [c3f,BkVec] = Cos3f(M,e,stop)
if nargin < 3 % default tolerance
    tol = 1e-14;
    nTol = tol;
    kMax = Inf;
elseif stop < 1 % stop is tolerance
    tol = stop;
    nTol = stop;
    kMax = Inf;
else            % stop is maximum iterations
    tol = 0;
    kMax = stop;
    nTol = 1e-14;
end
b = (1-sqrt(1-e^2))/e;
c3f = zeros(1,length(M));
B0 = 2*b^2*(b^5 -3*e*b^4 -2*b^3 +(2*e^3+12*e)*b^2 -(12*e^2+3)*b ...
        +2*e^3 +3*e)/(e^2*(b^2-1)^3); % constant term
Am = [1,-6*e,12*e^2+3,-8*e^3-12*e,12*e^2+3,-6*e,1];
m = [0:6].';

BkVec = [B0];

Bk = inf;
k = 1;
while abs(Bk) > tol && k < kMax % Iterate until Bk is small enough
    dBk = inf;
    n = 0;
    
    Bk = 1/4/(1-e^2)*besselj(n,-k*e)*Am*...
        (abs(m+n+k-2).*b.^abs(m+n+k-3) + e/sqrt(1-e^2)*b.^abs(m+n+k-2));
    n = 1;
    while abs(dBk) > nTol || n <= k + 3% Iterate until Bk converges
        dBk = 1/4/(1-e^2)*besselj(n,-k*e)*Am*...
        (abs(m+n+k-2).*b.^abs(m+n+k-3) + e/sqrt(1-e^2)*b.^abs(m+n+k-2))+...
        1/4/(1-e^2)*besselj(-n,-k*e)*Am*...
        (abs(m-n+k-2).*b.^abs(m-n+k-3) + e/sqrt(1-e^2)*b.^abs(m-n+k-2));
        
        Bk = Bk + dBk;
        n = n+1;
    end
    BkVec = [BkVec; Bk]; %#ok<AGROW>
    k = k+1;
end
for k = 0:length(BkVec)-1
    c3f = c3f + BkVec(k+1)*cos(k*M);
end
end

function [c4f,BkVec] = Cos4f(M,e,stop)
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
b = (1-sqrt(1-e^2))/e;
c4f = zeros(1,length(M));
B0 = b^3*(3*b^8-8*e*b^7-12*b^6+40*e*b^5+(8*e^4+24*e^2+18)*b^4-(96*e^3+152*e)*b^3+...
    (32*e^4+240*e^2+36)*b^2-(96*e^3+72*e)*b+8*e^4+24*e^2+3)/(e^3*(1-b^2)^5);% constant term
    
Am = [1,-8*e,24*e^2+4,-(32*e^3+24*e),(16*e^4+48*e^2+6),-(32*e^3+24*e),24*e^2+4,-8*e,1];
m = [0:8].'; %#ok<NBRAK>

BkVec = [B0];

Bk = inf;
k = 1;
while abs(Bk) > tol && k < kMax % Iterate until Bk is small enough
    dBk = inf;
    n = 0;
    
    Bk = 1/8/(1-e^2)^(3/2)*besselj(n,-k*e)*Am*...
        (abs(m+n+k-3).*abs(m+n+k-2)/2.*b.^abs(m+n+k-4) + ...
        3/2*e/sqrt(1-e^2)*abs(m+n+k-2).*b.^abs(m+n+k-3) + ...
        3/2*e^2/(1-e^2)*b.^abs(m+n+k-2));
    
    n = 1;
    while abs(dBk) > nTol || n <= k + 4% Iterate until Bk converges
        dBk = 1/8/(1-e^2)^(3/2)*besselj(n,-k*e)*Am*...
        (abs(m+n+k-3).*abs(m+n+k-2)/2.*b.^abs(m+n+k-4) + ...
        3/2*e/sqrt(1-e^2)*abs(m+n+k-2).*b.^abs(m+n+k-3) + ...
        3/2*e^2/(1-e^2)*b.^abs(m+n+k-2))...
        ...
        +1/8/(1-e^2)^(3/2)*besselj(-n,-k*e)*Am*...
        (abs(m-n+k-3).*abs(m-n+k-2)/2.*b.^abs(m-n+k-4) + ...
        3/2*e/sqrt(1-e^2)*abs(m-n+k-2).*b.^abs(m-n+k-3) + ...
        3/2*e^2/(1-e^2)*b.^abs(m-n+k-2));
        
        Bk = Bk + dBk;
        n = n+1;
    end
    BkVec = [BkVec; Bk];
    k = k+1;
end
for k = 0:length(BkVec)-1
    c4f = c4f + BkVec(k+1)*cos(k*M);
end
end

function [c5f,BkVec] = Cos5f(M,e,stop)
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
b = (1-sqrt(1-e^2))/e;
c5f = zeros(1,length(M));
B0 = -2*b^4*(2*b^11 -5*e*b^10 -12*b^9 +35*e*b^8 +28*b^7 ...
    +(8*e^5+40*e^3-90*e)*b^6 -(160*e^4+240*e^2+48)*b^5 ...
    +(72*e^5+760*e^3+510*e)*b^4 -(480*e^4+1120*e^2+110)*b^3 ...
    +(72*e^5+760*e^3+335*e)*b^2 -(160*e^4+240*e^2+20)*b +8*e^5 +40*e^3 +15*e)...
    /e^4/(1-b^2)^7;% constant term
    
Am = [1,-10*e,40*e^2+5,-(80*e^3+40*e),(80*e^4+120*e^2+10),...
    -(32*e^5+160*e^3+60*e),(80*e^4+120*e^2+10),-(80*e^3+40*e),40*e^2+5,...
    -10*e,1];
m = [0:10].';

BkVec = [B0];

Bk = inf;
k = 1;
while abs(Bk) > tol && k < kMax % Iterate until Bk is small enough
    dBk = inf;
    n = 0;
    
    Bk = 1/16/(1-e^2)^(2)*besselj(n,-k*e)*Am*...
        (abs(m+n+k-4).*abs(m+n+k-3).*abs(m+n+k-2)/6.*b.^abs(m+n+k-5) + ...
        e/sqrt(1-e^2)*abs(m+n+k-3).*abs(m+n+k-2).*b.^abs(m+n+k-4) + ...
        5/2*e^2/(1-e^2)*abs(m+n+k-2).*b.^abs(m+n+k-3) + ...
        5/2*e^3/(1-e^2)^(3/2)*b.^abs(m+n+k-2));
    
    n = 1;
    while abs(dBk) > nTol || n <= k + 5% Iterate until Bk converges
        dBk = 1/16/(1-e^2)^(2)*besselj(n,-k*e)*Am*...
        (abs(m+n+k-4).*abs(m+n+k-3).*abs(m+n+k-2)/6.*b.^abs(m+n+k-5) + ...
        e/sqrt(1-e^2)*abs(m+n+k-3).*abs(m+n+k-2).*b.^abs(m+n+k-4) + ...
        5/2*e^2/(1-e^2)*abs(m+n+k-2).*b.^abs(m+n+k-3) + ...
        5/2*e^3/(1-e^2)^(3/2)*b.^abs(m+n+k-2))...
        ...
        +1/16/(1-e^2)^(2)*besselj(-n,-k*e)*Am*...
        (abs(m-n+k-4).*abs(m-n+k-3).*abs(m-n+k-2)/6.*b.^abs(m-n+k-5) + ...
        e/sqrt(1-e^2)*abs(m-n+k-3).*abs(m-n+k-2).*b.^abs(m-n+k-4) + ...
        5/2*e^2/(1-e^2)*abs(m-n+k-2).*b.^abs(m-n+k-3) + ...
        5/2*e^3/(1-e^2)^(3/2)*b.^abs(m-n+k-2));
        
        Bk = Bk + dBk;
        n = n+1;
    end
    BkVec = [BkVec; Bk];
    k = k+1;
end
for k = 0:length(BkVec)-1
    c5f = c5f + BkVec(k+1)*cos(k*M);
end
end

function [c1s1f,BkVec] = CosfSinf(M,e,stop)
if nargin < 3 % default tolerance
    tol = 1e-14;
    nTol = tol;
    kMax = Inf;
elseif stop < 1 % stop is tolerance
    tol = stop;
    nTol = stop;
    kMax = Inf;
else            % stop is maximum iterations
    tol = 0;
    kMax = stop;
    nTol = 1e-14;
end
b = (1-sqrt(1-e^2))/e;
c1s1f = zeros(1,length(M));
Am = [-1,2*e,0,-2*e,1];
m = [0:4].';

BkVec = [];

Bk = inf;
k = 1;
while abs(Bk) > tol && k < kMax % Iterate until Bk is small enough
    dBk = inf;
    n = 0;
    
    Bk = -1/2*besselj(n,-k*e)*Am*...
        (b.^abs(m+n+k-2));
    n = 1;
    while abs(dBk) > nTol || n <= k+2% Iterate until Bk converges
        dBk = -1/2*besselj(n,-k*e)*Am*...
        (b.^abs(m+n+k-2))+...
        -1/2*besselj(-n,-k*e)*Am*...
        (b.^abs(m-n+k-2));
        
        Bk = Bk + dBk;
        n = n+1;
    end
    BkVec = [BkVec; Bk];
    k = k+1;
end
for k = 1:length(BkVec)
    c1s1f = c1s1f + BkVec(k)*sin(k*M);
end
end

function [c2s1f,BkVec] = Cos2fSinf(M,e,stop)
if nargin < 3 % default tolerance
    tol = 1e-14;
    nTol = tol;
    kMax = Inf;
elseif stop < 1 % stop is tolerance
    tol = stop;
    nTol = stop;
    kMax = Inf;
else            % stop is maximum iterations
    tol = 0;
    kMax = stop;
    nTol = 1e-14;
end
b = (1-sqrt(1-e^2))/e;
c2s1f = zeros(1,length(M));
Am = [-1,4*e,-4*e^2-1,0,4*e^2+1,-4*e,1];
m = [0:6].';

BkVec = [];

Bk = inf;
k = 1;
while abs(Bk) > tol && k < kMax % Iterate until Bk is small enough
    dBk = inf;
    n = 0;
    
    Bk = -1/4/sqrt(1-e^2)*besselj(n,-k*e)*Am*...
        (abs(m+n+k-2).*b.^abs(m+n+k-3) + ...
        e/sqrt(1-e^2)*b.^abs(m+n+k-2));
    n = 1;
    while abs(dBk) > nTol || n <= k+3% Iterate until Bk converges
        dBk = -1/4/sqrt(1-e^2)*besselj(n,-k*e)*Am*...
        (abs(m+n+k-2).*b.^abs(m+n+k-3) + ...
        e/sqrt(1-e^2)*b.^abs(m+n+k-2))+...
        -1/4/sqrt(1-e^2)*besselj(-n,-k*e)*Am*...
        (abs(m-n+k-2).*b.^abs(m-n+k-3) + ...
        e/sqrt(1-e^2)*b.^abs(m-n+k-2));
        
        Bk = Bk + dBk;
        n = n+1;
    end
    BkVec = [BkVec; Bk];
    k = k+1;
end
for k = 1:length(BkVec)
    c2s1f = c2s1f + BkVec(k)*sin(k*M);
end
end

function [c3s1f,BkVec] = Cos3fSinf(M,e,stop)
if nargin < 3 % default tolerance
    tol = 1e-14;
    nTol = tol;
    kMax = Inf;
elseif stop < 1 % stop is tolerance
    tol = stop;
    nTol = stop;
    kMax = Inf;
else            % stop is maximum iterations
    tol = 0;
    kMax = stop;
    nTol = 1e-14;
end
b = (1-sqrt(1-e^2))/e;
c3s1f = zeros(1,length(M));
Am = [-1,6*e,-12*e^2-2,8*e^3+6*e,0,-8*e^3-6*e,12*e^2+2,-6*e,1];
m = [0:8].';

BkVec = [];

Bk = inf;
k = 1;
while abs(Bk) > tol && k < kMax % Iterate until Bk is small enough
    dBk = inf;
    n = 0;
    
    Bk = -1/8/(1-e^2)*besselj(n,-k*e)*Am*...
        (abs(m+n+k-3).*abs(m+n+k-2)/2.*b.^abs(m+n+k-4) + ...
        3/2*e/sqrt(1-e^2)*abs(m+n+k-2).*b.^abs(m+n+k-3) + ...
        3/2*e^2/(1-e^2)*b.^abs(m+n+k-2));
    n = 1;
    while abs(dBk) > nTol || n <= k+4% Iterate until Bk converges
        dBk = -1/8/(1-e^2)*besselj(n,-k*e)*Am*...
        (abs(m+n+k-3).*abs(m+n+k-2)/2.*b.^abs(m+n+k-4) + ...
        3/2*e/sqrt(1-e^2)*abs(m+n+k-2).*b.^abs(m+n+k-3) + ...
        3/2*e^2/(1-e^2)*b.^abs(m+n+k-2))+ ...
        -1/8/(1-e^2)*besselj(-n,-k*e)*Am*...
        (abs(m-n+k-3).*abs(m-n+k-2)/2.*b.^abs(m-n+k-4) + ...
        3/2*e/sqrt(1-e^2)*abs(m-n+k-2).*b.^abs(m-n+k-3) + ...
        3/2*e^2/(1-e^2)*b.^abs(m-n+k-2));
        
        Bk = Bk + dBk;
        n = n+1;
    end
    BkVec = [BkVec; Bk];
    k = k+1;
end
for k = 1:length(BkVec)
    c3s1f = c3s1f + BkVec(k)*sin(k*M);
end
end

function [c4s1f,BkVec] = Cos4fSinf(M,e,stop)
if nargin < 3 % default tolerance
    tol = 1e-14;
    nTol = tol;
    kMax = Inf;
elseif stop < 1 % stop is tolerance
    tol = stop;
    nTol = stop;
    kMax = Inf;
else            % stop is maximum iterations
    tol = 0;
    kMax = stop;
    nTol = 1e-14;
end
b = (1-sqrt(1-e^2))/e;
c4s1f = zeros(1,length(M));
Am = [-1,8*e,-24*e^2-3,32*e^3+16*e,-16*e^4-24*e^2-2,0,16*e^4+24*e^2+2,...
    -32*e^3-16*e,24*e^2+3,-8*e,1];
m = [0:10].';

BkVec = [];

Bk = inf;
k = 1;
while abs(Bk) > tol && k < kMax % Iterate until Bk is small enough
    dBk = inf;
    n = 0;
    
    Bk = -1/16/(1-e^2)^(3/2)*besselj(n,-k*e)*Am*...
        (abs(m+n+k-4).*abs(m+n+k-3).*abs(m+n+k-2)/6.*b.^abs(m+n+k-5) + ...
        e/sqrt(1-e^2)*abs(m+n+k-3).*abs(m+n+k-2).*b.^abs(m+n+k-4) + ...
        5/2*e^2/(1-e^2)*abs(m+n+k-2).*b.^abs(m+n+k-3) + ...
        5/2*e^3/(1-e^2)^(3/2)*b.^abs(m+n+k-2));
    n = 1;
    while abs(dBk) > nTol || n <= k+5% Iterate until Bk converges
        dBk = -1/16/(1-e^2)^(3/2)*besselj(n,-k*e)*Am*...
        (abs(m+n+k-4).*abs(m+n+k-3).*abs(m+n+k-2)/6.*b.^abs(m+n+k-5) + ...
        e/sqrt(1-e^2)*abs(m+n+k-3).*abs(m+n+k-2).*b.^abs(m+n+k-4) + ...
        5/2*e^2/(1-e^2)*abs(m+n+k-2).*b.^abs(m+n+k-3) + ...
        5/2*e^3/(1-e^2)^(3/2)*b.^abs(m+n+k-2))+ ...
        -1/16/(1-e^2)^(3/2)*besselj(-n,-k*e)*Am*...
        (abs(m-n+k-4).*abs(m-n+k-3).*abs(m-n+k-2)/6.*b.^abs(m-n+k-5) + ...
        e/sqrt(1-e^2)*abs(m-n+k-3).*abs(m-n+k-2).*b.^abs(m-n+k-4) + ...
        5/2*e^2/(1-e^2)*abs(m-n+k-2).*b.^abs(m-n+k-3) + ...
        5/2*e^3/(1-e^2)^(3/2)*b.^abs(m-n+k-2));
        
        Bk = Bk + dBk;
        n = n+1;
    end
    BkVec = [BkVec; Bk];
    k = k+1;
end
for k = 1:length(BkVec)
    c4s1f = c4s1f + BkVec(k)*sin(k*M);
end
end

function f = fFourier(M,e,tol)
if nargin < 3
    tol = 1e-14;
end
b = (1-sqrt(1-e^2))/e;
f = nan(1,length(M));
for iM = 1:length(M)
    f(iM) = M(iM);
    df = inf;
    k = 1;
    while abs(df) > tol
        dBk = inf;
        n = 1;
        Bk = 2/k*besselj(0,-k*e)*b^abs(k);
        while abs(dBk) > tol
            dBk = 2/k*(besselj(n,-k*e)*b^abs(n+k) + ...
                besselj(-n,-k*e)*b^abs(-n+k));
            
            Bk = Bk + dBk;
            n = n+1;
        end
        df = Bk*sin(k*M(iM));
        f(iM) = f(iM) + df;
        k = k+1;
    end
end
end