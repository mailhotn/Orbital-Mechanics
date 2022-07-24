clear
%% Scale & Tolerance
primary = earth();
% True Scale
mu = primary.mu;
Re = primary.Re;
J2 = primary.J2;
dScale = 1;
tScale = 1;
% Normalized
% dScale = sqrt(primary.J2)*primary.Re;
% tScale = sqrt(dScale^3/primary.mu);
% mu = 1;
% Re = primary.Re/dScale;
% J2 = 1/Re^2;

imagTol = 1e-8;
%% Initial Conditions
sma = 10000/dScale;
ecc = 0.2;
inc = 60;
ran = 30;
aop = 30;
man = 0;
f = me2ta(man,ecc);
nT = 1000;
oeC = nan(6,nT);
oeW = nan(6,nT);


%% Coordinate Switch
radQ = sma*((1-ecc^2)/(1+ecc*cosd(f)));     % r
aolQ = pi/180*(aop + f);                  % theta
ranQ = pi/180*ran;                        % nu
vraP = sqrt(mu/sma/(1-ecc^2))*ecc*sind(f);% R
amoP = sqrt((1-ecc^2)*sma*mu);            % Theta
amzP = amoP*cosd(inc);                    % N

%% Solution
X = mu*J2*Re^2*(0.5-1.5*amzP^2/amoP^2);

h = 0.5*vraP^2 + 0.5*amoP^2/radQ^2 - mu/radQ + ...
    0.25*mu*J2*Re^2/radQ^3*(1-3*amzP^2/amoP^2); % initial energy
% hUp = (amoP^6+9*mu*X*amoP^2+(amoP^4+6*mu*X)^(3/2))/(27*X^2);
% hLow = (amoP^6+9*mu*X*amoP^2-(amoP^4+6*mu*X)^(3/2))/(27*X^2);
% iS = 1:3;
% l = acos((27*h/X-9*mu*amoP^2/X^2-amoP^6/X^3)/(amoP^4/X^2+6*mu/X)^(3/2));
% l2 = acos((-27*h/X+9*mu*amoP^2/X^2+amoP^6/X^3)/(amoP^4/X^2+6*mu/X)^(3/2));
% sSol = -1/3*amoP^2/X + 2/3*sqrt(amoP^4/X^2 + 6*mu/X)*cos((l+2*pi*iS)/3);
% sSol2 = -1/3*amoP^2/X - 2/3*sqrt(amoP^4/X^2 + 6*mu/X)*cos((l2+2*pi*(iS-2))/3);
% Numerical root finding - most accurate
sSol3 = roots([1,amoP^2/X,-2*mu/X,-2*h/X]);
sSol = sort(sSol3);
s1 = sSol(1);
s2 = sSol(2);
s3 = sSol(3);
% Diagnostics - energy at apses
h1 = 0.5*amoP^2*s1^2 - mu*s1 + 0.25*mu*J2*Re^2*s1^3*(1-3*amzP^2/amoP^2);
h2 = 0.5*amoP^2*s2^2 - mu*s2 + 0.25*mu*J2*Re^2*s2^3*(1-3*amzP^2/amoP^2);
h3 = 0.5*amoP^2*s3^2 - mu*s3 + 0.25*mu*J2*Re^2*s3^3*(1-3*amzP^2/amoP^2);
if X < 0    
    k0 = (s2-s1)/(s3-s1);
    k1 = 1/k0;
    n0 = (s2-s1)/s1;
    n1 = n0*k1;
    sV = linspace(1/radQ,s1,nT);
    z0 = (sV-s1)/(s2-s1);
    if max(z0-1) < imagTol % remove small imaginary stuff
        z0(z0>1)=1;
    else
        error('Apsis error too large!')
    end
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
    if max(z0-1) < imagTol % remove small imaginary stuff
        z0(z0>1)=1;
    else
        error('Apsis error too large!')
    end
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
hVec = 0.5*oeW(4,:).^2 + 0.5*oeW(5,:).^2.*sV.^2 - mu.*sV + ...
    0.25*mu*J2*Re^2.*sV.^3.*(1-3*oeW(6,:).^2./oeW(5,:).^2);
if max(abs(hVec-h)) < imagTol
    oeW = real(oeW);
else
    error('Energy error too large!')
end
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

figure(7)
plot(t,fVec)