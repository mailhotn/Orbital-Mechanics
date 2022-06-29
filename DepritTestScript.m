clear
primary = earth();
mu = primary.mu;
Re = primary.Re;
J2 = primary.J2;
sma = 10000;
ecc = 0.2;
inc = 40;
ran = 30;
aop = 30;
man = 1;
f = me2ta(man,ecc);
nT = 1000;
oeC = nan(6,nT);
oeW = nan(6,nT);
%% Coordinate Switch
radQ = sma*(1-ecc^2)/(1+ecc*cosd(f));     % r
aolQ = pi/180*(aop + f);                  % theta
ranQ = pi/180*ran;                        % nu
vraP = sqrt(mu/sma/(1-ecc^2))*ecc*sind(f);% R
amoP = sqrt((1-ecc^2)*sma*mu);            % Theta
amzP = amoP*cosd(inc);                    % N

%% Solution
X = mu*J2*Re^2*(0.5-1.5*amzP^2/amoP^2);
iS = 1:3;
h = 0.5*vraP^2 + 0.5*amoP^2/radQ^2 - mu/radQ + ...
    0.25*mu*J2*Re^2/radQ^3*(1-3*amzP^2/amoP^2);
hUp = (amoP^6+9*mu*X*amoP^2+(amoP^4+6*mu*X)^(3/2))/(27*X^2);
hLow = (amoP^6+9*mu*X*amoP^2-(amoP^4+6*mu*X)^(3/2))/(27*X^2);
l = acos((27*h/X-9*mu*amoP^2/X^2-amoP^6/X^3)/(amoP^4/X^2+6*mu/X)^(3/2));
sSol = -1/3*amoP^2/X + 2/3*sqrt(amoP^4/X^2 + 6*mu/X)*cos((l+2*pi*iS)/3);
s1 = sSol(1);
s2 = sSol(2);
s3 = sSol(3);
if X < 0    
    k0 = (s2-s1)/(s3-s1);
    k1 = 1/k0;
    n0 = (s2-s1)/s1;
    n1 = n0*k1;
    sV = linspace(1/radQ,s1,nT);
    z0 = (sV-s1)/(s2-s1);
    z1 = z0*k0;
    phi = asin(sqrt(z0));
    
    % RAAN solution
    Iv0 = 2*sqrt(s3/-X)*(sqrt(s3/(s3-s1))*ellipticK(k0)-...
        sqrt((s3-s1)/s3)*ellipticE(k0));
    Iv = Iv0 - 2*sqrt(s3/-X)*(sqrt(s3/(s3-s1))*ellipticF(phi,k0)-...
        sqrt((s3-s1)/s3)*ellipticE(phi,k0));
    v0 = -3/2*mu*J2*Re^2*amzP/amoP^2*Iv0;
    
    % AOL solution
    Ith0 = 2/(sqrt(-X)*sqrt(s3-s1))*ellipticK(k0);
    Ith = Ith0 - 2/(sqrt(-X)*sqrt(s3-s1))*ellipticF(phi,k0);
    th0 = amoP*Ith0 - v0*amzP/amoP;
    
    % Time Solution
    T0 = 1/(sqrt(-X)*s1^2*sqrt(s2-s1))*sqrt(k0)/(1+n0)*...
        (((3+2*n0)*k0+(2+n0)*n0)/(k0+n0)*ellipticPi(-n0,k0) + ...
        n0/(k0+n0)*ellipticE(k0) - ellipticK(k0));
    It = T0 - 1/(sqrt(-X)*s1^2*sqrt(s2-s1))*sqrt(k0)/(1+n0)*...
        (n0/2*n0/(k0+n0)*sqrt(1-k0*sin(phi).^2)./(1+n0*sin(phi).^2).*sin(2*phi)...
        +((3+2*n0)*k0+(2+n0)*n0)/(k0+n0)*ellipticPi(-n0,phi,k0) ...
        + n0/(k0+n0)*ellipticE(phi,k0) -ellipticF(phi,k0));
    t = It;
    
elseif X > 0
    k0 = 1-(s2-s1)/(s3-s1);
    k1 = 1/k0;
    n0 = (s3-s2)/s3;
    n1 = n0*k1;
    sV = linspace(1/radQ,s2,nT);
    z0 = (s3-sV)/(s3-s2);
    z1 = z0*k0;
    phi = asin(sqrt(z0));
    
    % RAAN solution
    Iv0 = 2*sqrt(s1/X)*(sqrt(s1/(s3-s1))*ellipticK(k0)-...
        sqrt((s3-s1)/s1)*ellipticE(k0));
    Iv = 2*sqrt(s1/X)*(sqrt(s1/(s3-s1))*ellipticF(phi,k0)-...
        sqrt((s3-s1)/s1)*ellipticE(phi,k0));
    v0 = -3/2*mu*J2*Re^2*amzP/amoP^2*Iv0;
    
    % AOL solution
    Ith0 = 2/(sqrt(X)*sqrt(s3-s1))*ellipticK(k0);
    Ith = 2/(sqrt(X)*sqrt(s3-s1))*ellipticF(phi,k0);
    th0 = amoP*Ith0 - v0*amzP/amoP;
    
    % Time Solution
    T0 = 1/(sqrt(X)*s3^2*sqrt(s3-s2))*sqrt(k0)/(1-n0)*...
        (((3-2*n0)*k0-(2*k0-n0)*n0)/(k0-n0)*ellipticPi(n0,k0) - ...
        n0/(k0-n0)*ellipticE(k0) - ellipticK(k0));
    It = 1/(sqrt(X)*s3^2*sqrt(s3-s2))*sqrt(k0)/(1-n0)*...
        (n0/2*n0/(k0-n0)*sqrt(1-k0*sin(phi).^2)./(1-n0*sin(phi).^2).*sin(2*phi)...
        +(3*k0-2*n0-(2*k0-n0)*n0)/(k0-n0)*ellipticPi(n0,phi,k0) ...
        - n0/(k0-n0)*ellipticE(phi,k0) -ellipticF(phi,k0));
    t = It;
    
else
end

% Finish solution
oeW(1,:) = 1./sV; % radial position
oeW(2,:) = aolQ + 1.5*mu*J2*Re^2*amzP^2/amoP^3*Iv + amoP*Ith;
oeW(3,:) = ranQ -1.5*mu*J2*Re^2*amzP/amoP^2*Iv;
oeW(4,:) = sqrt(2*h + 2*mu.*sV - amoP^2.*sV.^2 - X.*sV.^3);
oeW(5,:) = amoP;
oeW(6,:) = amzP;
%% Convert to Conventional
oeC(1,:) = -mu*oeW(1,:).^2/(oeW(1,:).^2.*oeW(4,:).^2+...
    oeW(5,:).^2-2*mu*oeW(1,:));
oeC(2,:) = sqrt(1-oeW(5,:).^2./(mu*oeC(1,:)));
oeC(3,:) = acosd(oeW(6,:)./oeW(5,:));
oeC(4,:) = 180/pi*oeW(3,:);
fVec = atan2(oeW(5,:).*oeW(4,:)./(mu*oeC(2,:)),...
    (oeW(5,:).^2-mu*oeW(1,:))./(mu*oeW(1,:).*oeC(2,:)));
oeC(5,:) = 180/pi*(oeW(2,:) - fVec);
oeC(6,:) = ta2me(fVec*180/pi,oeC(2,:));



%% Plot
figure(1)
plot(t,oeW(1,:))

figure(2)
plot(t,oeW(2,:))

figure(3)
plot(t,oeW(3,:))

figure(4)
plot(t,oeW(4,:))

figure(5)
plot(t,oeC(5,:))

figure(6)
plot(t,oeW(6,:))

