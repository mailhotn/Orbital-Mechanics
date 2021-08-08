e = 0.1;
M = 0:0.001:2*pi;
tol = 1e-14;
%% Calculate N-R, Battin

fFour = fFourier(M,e,tol);

tic
fTrue = pi/180*me2ta(180/pi*M,e,tol);
c2fTrue = cos(fTrue).^2;
c3fTrue = cos(fTrue).^3;
c4fTrue = cos(fTrue).^4;
toc

tic
[c2fFourier, BkVec2] = Cos2f(M,e,tol);
toc
tic
[c3fFourier, BkVec3] = Cos3f(M,e,tol);
toc
tic
[c4fFourier, BkVec4] = Cos4f(M,e,tol);
toc


errorc2f = norm(c2fTrue-c2fFourier)/length(c2fTrue)
errorc3f = norm(c3fTrue-c3fFourier)/length(c3fTrue)
errorc4f = norm(c4fTrue-c4fFourier)/length(c4fTrue)
%% Plot Stuff
figure(1) % M & E
plot(M,fTrue,M,fFour,'--','linewidth',2)
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
yticks(0:pi/2:2*pi)
yticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
axis([0 2*pi 0 2*pi])

figure(2)
plot(M,c2fTrue,M,c2fFourier,'--','linewidth',2)
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('M')
axis([0 2*pi -1 1])
legend('N-R','Fourier','location','best')

figure(3)
plot(M,c3fTrue,M,c3fFourier,'--','linewidth',2)
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('M')
axis([0 2*pi -1 1])
legend('N-R','Fourier','location','best')

figure(4)
plot(M,c4fTrue,M,c4fFourier,'--','linewidth',2)
xticks(0:pi/2:2*pi)
xticklabels({'0','\pi/2','\pi','3\pi/2','2\pi'})
xlabel('M')
axis([0 2*pi -1 1])
legend('N-R','Fourier','location','best')



%%
function [c2f,BkVec] = Cos2f(M,e,stop)
if nargin < 3 % default tolerance
    tol = 1e-14;
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
c2f = nan(1,length(M));
B0 = b*(2*e^2-4*e*b+b^2+1)/e/(1-b^2); % constant term
BkVec = [];

Bk = inf;
k = 1;
while abs(Bk) > tol && k < kMax% Iterate until Bk is small enough

    dBk = inf;
    n = 1;
    Bk = 1/2/sqrt(1-e^2)*besselj(0,-k*e)*(b^abs(k+2)-4*e*b^abs(k+1)+...
        (4*e^2+2)*b^abs(k)-4*e*b^abs(k-1)+b^abs(k-2));
    while abs(dBk) > nTol || n < k% Iterate until Bk converges
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
c2f(:) = B0;
for k = 1:length(BkVec)
    c2f = c2f + BkVec(k)*cos(k*M);
end
end

function [c3f,BkVec] = Cos3f(M,e,stop)
if nargin < 3 % default tolerance
    tol = 1e-14;
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
c3f = nan(1,length(M));
B0 = 2*b^2*(b^5 -3*e*b^4 -2*b^3 +(2*e^3+12*e)*b^2 -(12*e^2+3)*b ...
        +2*e^3 +3*e)/(e^2*(b^2-1)^3); % constant term
Am = [1,-6*e,12*e^2+3,-8*e^3-12*e,12*e^2+3,-6*e,1];
m = [0:6].';

BkVec = [];

Bk = inf;
k = 1;
while abs(Bk) > tol && k < kMax % Iterate until Bk is small enough
    dBk = inf;
    n = 0;
    
    Bk = 1/4/(1-e^2)*besselj(n,-k*e)*Am*...
        (abs(m+n+k-2).*b.^abs(m+n+k-3) + e/sqrt(1-e^2)*b.^abs(m+n+k-2));
    n = 1;
    while abs(dBk) > nTol || n < k% Iterate until Bk converges
        dBk = 1/4/(1-e^2)*besselj(n,-k*e)*Am*...
        (abs(m+n+k-2).*b.^abs(m+n+k-3) + e/sqrt(1-e^2)*b.^abs(m+n+k-2))+...
        1/4/(1-e^2)*besselj(-n,-k*e)*Am*...
        (abs(m-n+k-2).*b.^abs(m-n+k-3) + e/sqrt(1-e^2)*b.^abs(m-n+k-2));
        
        Bk = Bk + dBk;
        n = n+1;
    end
    BkVec = [BkVec; Bk];
    k = k+1;
end
c3f(:) = B0;
for k = 1:length(BkVec)
    c3f = c3f + BkVec(k)*cos(k*M);
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
c4f = nan(1,length(M));
B0 = b^3*(3*b^8-8*e*b^7-12*b^6+40*e*b^5+(8*e^4+24*e^2+18)*b^4-(96*e^3+152*e)*b^3+...
    (32*e^4+240*e^2+36)*b^2-(96*e^3+72*e)*b+8*e^4+24*e^2+3)/(e^3*(1-b^2)^5);% constant term
    
Am = [1,-8*e,24*e^2+4,-(32*e^3+24*e),(16*e^4+48*e^2+6),-(32*e^3+24*e),24*e^2+4,-8*e,1];
m = [0:8].';

BkVec = [];

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
    while abs(dBk) > nTol || n < k % Iterate until Bk converges
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
c4f(:) = B0;
for k = 1:length(BkVec)
    c4f = c4f + BkVec(k)*cos(k*M);
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