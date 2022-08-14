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

imagTol = 1e-12;
nOrb = 1; % number of orbits
nTArc = 40; % number of time-steps per half orbit
nTime = nTArc*nOrb*2; % number of time steps
%% Initial Conditions
sma = 25000/dScale;
ecc = 0.5;
inc = 1;
% inc = 180-acosd(1/sqrt(3)); % singularity test
ran = 50;
aop = 359;
man = 0;

oeC = nan(nTime,6);
oeW = nan(nTime,6);

%% Bad IC
dataStruct = open('DepritTErrorIC.mat');
oeErr = dataStruct.oeErr;
iErr = 1;
sma = oeErr(iErr,1);
ecc = oeErr(iErr,2);
inc = oeErr(iErr,3);
ran = oeErr(iErr,4);
aop = oeErr(iErr,5);
man = oeErr(iErr,6);
%% Coordinate Switch
f = me2ta(man,ecc);
radQ = sma*((1-ecc^2)/(1+ecc*cosd(f)));   % r
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
p = [1,amoP^2/X,-2*mu/X,-2*h/X];
A = diag(ones(2,1),-1);
A(1,:) = -p(2:4)./p(1);

if cond(A) < 1/imagTol % Check cond number of companion matrix
    % Well conditioned, i is not too close to critical - numerical errors
    % will be reasonable
    sSol3 = roots(p);
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
        if f < 180 % Ascending IC
            tS = pi*(s2-1/radQ)/(s2-s1)+linspace(0,2*pi*nOrb,nTime).';
%             sVec = s1+(s2-s1)/2*(1-sawtooth(tS,0.5));
            sVec = s1+(s2-s1)/2*(1+cos(tS));
        else % Descending IC
            tS = pi*(1/radQ-s1)/(s2-s1)+linspace(0,2*pi*nOrb,nTime).';
%             sVec = s1+(s2-s1)/2*(1+sawtooth(tS,0.5));
            sVec = s1+(s2-s1)/2*(1-cos(tS));
        end
        signR = square(tS);
        
        z0 = (sVec-s1)/(s2-s1);
        if max(z0-1) < imagTol && min(z0) > -imagTol% remove small imaginary stuff
            z0(z0>1) = 1;
            z0(z0<0) = 0;
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
        if f < 180 % Ascending IC
            tS = pi*(s3-1/radQ)/(s3-s2)+linspace(0,2*pi*nOrb,nTime).';
%             sVec = s2+(s3-s2)/2*(1-sawtooth(tS,0.5));
            sVec = s2+(s3-s2)/2*(1+cos(tS));
        else % Descending IC
            tS = pi*(1/radQ-s2)/(s3-s2)+linspace(0,2*pi*nOrb,nTime).';
%             sVec = s2+(s3-s2)/2*(1+sawtooth(tS,0.5));
            sVec = s2+(s3-s2)/2*(1-cos(tS));
        end
        signR = square(tS);
        z0 = (s3-sVec)/(s3-s2);
        if max(z0-1) < imagTol && min(z0) > -imagTol% remove small imaginary stuff
            z0(z0>1) = 1;
            z0(z0<0) = 0;
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
        v0 = -3/2*mu*J2*Re^2*amzP/amoP^2*Iv0; %v*
        
        % AOL solution
        Ith0 = 2/(sqrt(X)*sqrt(s3-s1))*ellipticK(k0);
        Ith = 2/(sqrt(X)*sqrt(s3-s1))*ellipticF(phi,k0);
        th0 = amoP*Ith0 - v0*amzP/amoP; %th*
        
        % Time Solution
        T0 = 1/(sqrt(X)*s3^2*sqrt(s3-s2))*sqrt(k0)/(1-n0)*...
            (((3-2*n0)*k0-(2-n0)*n0)/(k0-n0)*ellipticPi(n0,k0) - ...
            n0/(k0-n0)*ellipticE(k0) - ellipticK(k0));
        It = 1/(sqrt(X)*s3^2*sqrt(s3-s2))*sqrt(k0)/(1-n0)*...
            (n0/2*n0/(k0-n0)*sqrt(1-k0*sin(phi).^2)./(1-n0*sin(phi).^2).*sin(2*phi)...
            +(3*k0-2*n0-(2*k0-n0)*n0)/(k0-n0)*ellipticPi(n0,phi,k0) ...
            - n0/(k0-n0)*ellipticE(phi,k0) -ellipticF(phi,k0));
        t = It;
    end
    
    % unwrap Time - T0 is half period
    t = T0/(4*pi)*unwrap(4*pi/T0*(t).*signR);
    t = t-t(1);
    % Finish solution
    rVec = 1./sVec; % radial position
    % Fix Velocity sign
    RVec = signR.*sqrt(2*h + 2*mu.*sVec - amoP^2.*sVec.^2 - X.*sVec.^3);
    hVec = 0.5*RVec.^2 + 0.5*amoP.^2.*sVec.^2 - mu.*sVec + ...
        0.25*mu*J2*Re^2.*sVec.^3.*(1-3*amzP.^2./amoP.^2);
    if max(abs(hVec-h)) < imagTol
        RVec = real(RVec);
    else
        error('Energy error too large!')
        
    end
    fVec = unwrap(atan2(amoP.*RVec./(mu*ecc),...
    (amoP.^2-mu*rVec)./(mu*rVec.*ecc)));
else
    % Companion matrix is ill-conditioned, i is near critical,
    % and numerical errors would be large.
    % Keplerian Orbit
    T0 = 2*pi*sqrt(sma^3/mu)/2;
    t = linspace(0,T0*2*nOrb,nTime);
    manVec = man + 180/pi*sqrt(mu/sma^3)*t; % deg
    fVec = unwrap(pi/180*me2ta(manVec,ecc)).'; % rad
    % RAAN solution
    Iv = mu/amoP^3*(fVec+ecc*sin(fVec));
    v0 = -1.5*J2*mu^2*Re^2*amzP/amoP^5*pi;
    % AOL solution
    Ith = fVec/amoP;
    th0 = pi*(1+1.5*J2*mu^2*Re^2*amzP^2/amoP^6);
    % Other stuff
    rVec = sma*(1-ecc^2)./(1+ecc*cos(fVec));
    sVec = 1./rVec;
    RVec = sqrt(mu/sma/(1-ecc^2))*ecc*sin(fVec);
    signR = 1;
    vVec = sqrt(mu/sma/(1-ecc^2))*sqrt(1+ecc^2+2*ecc*cos(fVec));
    
    hVec = 0.5*vVec.^2 - mu./rVec +...
        0.25*mu*J2*Re^2./rVec.^3*(1-3*amzP^2/amoP^3); % should be 0
end
% Finish Solution
oeW(:,1) = rVec;
% unwrap AOL
dAol = 1.5*mu*J2*Re^2*amzP^2/amoP^3*Iv + amoP*Ith;
oeW(:,2) = th0/(2*pi)*unwrap(2*pi/th0*dAol.*signR)+aolQ;
% unwrap RAAN
dRan = -1.5*mu*J2*Re^2*amzP/amoP^2*Iv;
oeW(:,3) = v0/(2*pi)*unwrap(2*pi/v0*dRan.*signR)+ranQ;
oeW(:,4) = RVec;
oeW(:,5) = amoP;
oeW(:,6) = amzP;


%% Convert to Conventional
oeC(:,1) = -dScale*mu*oeW(:,1).^2./(oeW(:,1).^2.*oeW(:,4).^2+...
    oeW(:,5).^2-2*mu*oeW(:,1));
oeC(:,2) = sqrt(1-oeW(:,5).^2./(mu*oeC(:,1)/dScale));
oeC(:,3) = acosd(oeW(:,6)./oeW(:,5));
oeC(:,4) = 180/pi*oeW(:,3);
oeC(:,5) = 180/pi*(oeW(:,2) - fVec);
oeC(:,6) = ta2me(fVec*180/pi,oeC(:,2));
t = t*tScale;


%% Numerical Integration
oeI = [sma*dScale,ecc,inc,ran,aop,man];
Sat = SingleSat(oeI);
Prop = Propagator(Sat);
T = 2*pi*sqrt(sma^3/mu)*tScale;

[t2,oe] = Prop.PropOeOsc(t);


%% Plot
figure(1)
plot(t/T0/2,oeC(:,1))
xlabel('Orbits')
xlim([0,nOrb])
hold on
plot(t2/T,oe(:,1),'--')
hold off


figure(2)
plot(t/T0/2,oeC(:,2))
xlabel('Orbits')
xlim([0,nOrb])
hold on
plot(t2/T,oe(:,2),'--')
hold off


figure(3)
plot(t/T0/2,oeC(:,3))
xlabel('Orbits')
xlim([0,nOrb])
hold on
plot(t2/T,oe(:,3),'--')
hold off


figure(4)
plot(t/T0/2,oeC(:,4))
xlabel('Orbits')
xlim([0,nOrb])
hold on
plot(t2/T,oe(:,4),'--')
hold off


figure(5)
plot(t/T0/2,oeC(:,5))
xlabel('Orbits')
xlim([0,nOrb])
hold on
plot(t2/T,oe(:,5),'--')
hold off


figure(6)
plot(t/T0/2,oeC(:,6))
xlabel('Orbits')
xlim([0,nOrb])
hold on
plot(t2/T,oe(:,6),'--')
hold off


figure(7)
plot(t)
xlabel('Orbits')
% xlim([0,nOrb])




