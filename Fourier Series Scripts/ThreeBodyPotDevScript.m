e = 0.1;
a = 7000;

M = linspace(0,2*pi,10000);
f = me2ta(M*180/pi,e)*pi/180;
r = (1-e^2)*a./(1+e*cos(f));
r2c2fNum = r.^2.*cos(f).^2;
r2csfNum = r.^2.*cos(f).*sin(f);
kMax = 10;

aVec = [e, -4*e^2-2, 4*e^3+11*e, -16*e^2-4, 4*e^3+11*e, -4*e^2-2, e];
bVec = [e, -2*e^2-2, 5*e, 0, -5*e, 2*e^2+2, -e];
AVec = nan(1,kMax);
BVec = nan(1,kMax);
r2c2f = a^2*(4*e^2+1)/2*ones(size(M));
r2csf = zeros(size(M));
for k = 1:kMax
    nMin = -k-3;
    nMax = 3-k;
    Ak = 0;
    Bk = 0;
    for n = nMin:nMax
        m = 3-k-n;
        Ak = Ak + aVec(m+1)*besselj(n,-k*e);
        Bk = Bk + bVec(m+1)*besselj(n,-k*e);
    end
    Ak = -Ak*a^2/4;
    Bk = -Bk*a^2*sqrt(1-e^2)/4;
    AVec(k) = Ak;
    r2c2f = r2c2f + Ak*cos(k*M);
    r2csf = r2csf + Bk*sin(k*M);
end

figure(1)
plot(M,(r2c2fNum-r2c2f)./r2c2fNum)
figure(2)
plot(M,r2csfNum,M,r2csf)