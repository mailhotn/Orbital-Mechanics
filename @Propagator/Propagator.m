classdef Propagator < handle &  matlab.mixin.CustomDisplay
    % Propagator is a class that defines an orbit propagator for use with
    % a constellation object.
    properties(SetAccess = protected) % propagation
        relTol
        absTol
        Con     % Constellation object to propagate
    end
    
    methods
        function P = Propagator(Constellation, relTol, absTol)
            switch nargin
                case 0 % default tolerances, GPS constellation
                    relTol     = 1e-8;
                    absTol     = 1e-9;
                    Constellation = WalkerConstellation;
                case 1 % default tolerances
                    relTol     = 1e-8;
                    absTol     = 1e-9;
                case 3 % arbitrary orbit & tolerances
                    
                otherwise
                    error('Wrong number of input arguments')
            end
            P.relTol = relTol;
            P.absTol = absTol;
            P.Con    = Constellation;
        end
        
        function [Time, X] = PropEciTb(P,T)
            % Propagate for time T in ECI frame with no perturbations
            % This should be replaced with Analytical Solution
            % Complexity ~O(T)
            opts = odeset('reltol',P.relTol,'abstol',P.absTol);
            IC = reshape(P.Con.InitialStateEci,[6*P.Con.nSats,1]);
            [Time, X] = ode45(@P.DynEciTb,T,IC,opts);
        end
        
        function [Time, X] = PropEciJ2(P,T)
            % Propagate for time T in ECI frame with J2 perturbation
            % Complexity ~O(T)
            opts = odeset('reltol',P.relTol,'abstol',P.absTol);
            IC = reshape(P.Con.InitialStateEci,[6*P.Con.nSats,1]);
            [Time, X] = ode45(@P.DynEciJ2,T,IC,opts);
        end
        
        function [Time, X] = PropOeMean(P,T)
            % Basically useless, propagates equations of motion with linear
            % solution
            opts = odeset('reltol',P.relTol,'abstol',P.absTol);
            IC = reshape(P.Con.InitialOeMean,[6*P.Con.nSats,1]);
            [Time, X] = ode45(@P.DynOeMean,T,IC,opts);
        end
        
        function [Time, X] = PropOeMeanFast(P,T)
            % Directly calculates linear solution without numerical
            % integration
            IC = reshape(P.Con.InitialOeMean,[6*P.Con.nSats,1]);
            % move stuff around
            order = 6*P.Con.nSats;
            X2 = reshape(IC,[6,P.Con.nSats]);
            aV = X2(1,:);
            eV = X2(2,:);
            iV = X2(3,:);
            a = reshape(repmat(aV,6,1),order,1);
            e = reshape(repmat(eV,6,1),order,1);
            i = reshape(repmat(iV,6,1),order,1);
            % Derived values
            p = a.*(1-e.^2);
            n = sqrt(P.Con.mu./a.^3);
            eta = sqrt(1-e.^2);
            % Eq of motion
            dO = 180/pi*repmat([0,0,0,1,0,0].',P.Con.nSats,1).*...
                (-3/2*P.Con.J2.*(P.Con.Re./p).^2.*n.*cosd(i));
            dw = 180/pi*repmat([0,0,0,0,1,0].',P.Con.nSats,1).*...
                (3/4*P.Con.J2.*(P.Con.Re./p).^2.*n.*(5*cosd(i).^2-1));
            dM = 180/pi*repmat([0,0,0,0,0,1].',P.Con.nSats,1).*...
                (3/4*P.Con.J2.*(P.Con.Re./p).^2.*n.*eta.*(3*cosd(i).^2-1) + n);
            dX = dO + dw + dM;
            % "Propagate"
            X = zeros(order,length(T));
            X(:,1) = IC;
            X(:,2:end) = X(:,1) + dX*T(2:end);
            X = X.';
            Time = T;
        end
        
        function [Time, X] = PropOeOsc(P,T)
            % Use this I think
            opts = odeset('reltol',P.relTol,'abstol',P.absTol);
            OE = P.Con.InitialOeOsc;
            OE(3:end,:) = OE(3:end,:)*pi/180;
            IC = reshape(OE,[6*P.Con.nSats,1]);
            [Time, X] = ode45(@P.DynOeOsc,T,IC,opts);
            inddeg = reshape((3:6).'+(0:P.Con.nSats-1)*6,4*P.Con.nSats,1);
            indWrap = reshape((3:5).'+(0:P.Con.nSats-1)*6,3*P.Con.nSats,1);
            X(:,inddeg) = (180/pi*X(:,inddeg));
            X(:,indWrap) = wrapTo360(X(:,indWrap));
        end
        
        function [Time, X] = PropOeOsc2(P,T)
            % probably worse?
            opts = odeset('reltol',P.relTol,'abstol',P.absTol);
            OE = P.Con.InitialOeOsc;
            OE(3:end,:) = OE(3:end,:)*pi/180;
            IC = reshape(OE,[6*P.Con.nSats,1]);
            [Time, X] = ode45(@P.DynOeOsc2,T,IC,opts);
            inddeg = reshape((3:6).'+(0:P.Con.nSats-1)*6,4*P.Con.nSats,1);
            X(:,inddeg) = wrapTo360(180/pi*X(:,inddeg));
        end
        
        function [Time, X] = PropEciJ3(P,T)
            % Propagate for time T in ECI frame with J2 & J3 perturbations
            % Complexity ~O(T)
            opts = odeset('reltol',P.relTol,'abstol',P.absTol);
            IC = reshape(P.Con.InitialStateEci,[6*P.Con.nSats,1]);
            [Time, X] = ode45(@P.DynEciJ3,T,IC,opts);
        end
        
        function [Time, X] = PropOeFourierNum(P,T,kMax)
            % Propagate for time T using Fourier series of LPE
            % Very stupid
            % Currently only works for one satellite
            opts = odeset('reltol',P.relTol,'abstol',P.absTol);
            IC = reshape(P.Con.InitialOeOsc,[6*P.Con.nSats,1]);
            IC(3:end) = IC(3:end)*pi/180;
            [Time, X] = ode45(@(T,X) P.DynOeFourier(T,X,kMax),T,IC,opts);
            X(:,3:end) = X(:,3:end)*180/pi;
            X = X.';
            X(3:end,:) = wrapTo360(X(3:end,:));
        end
        
        function [Time, X] = PropOeFourier2(P,T,kMax)
            % Refactoring Bullshit
            [Time, X] = PropOeFourier(P,T,kMax);
        end
        
        function [Time, X] = PropOeFourier(P,T,kMax)
            % Propagate for time T using Fourier LPE
            % Assume constant coefficients, assume M is not affected by J2
            
            % Handle Initial conditions
            IC = reshape(P.Con.InitialOeOsc,[6*P.Con.nSats,1]);
            IC(3:end) = IC(3:end)*pi/180;
            a = IC(1);
            
            % Get mean a - nonsingular, somewhat weird
            e = IC(2);
            i = IC(3);
            aop = IC(5);
            eta = sqrt(1-e^2);
            f = pi/180*me2ta(IC(6)*180/pi,e);
            a_r = (1+e.*cos(f))./eta.^2;
            g2 = -P.Con.primary.J2/2*(P.Con.primary.Re/a)^2;
            a = a + a*g2*((3*cos(i)^2-1).*(a_r^3 - 1/eta^3) ...
                + 3*(1-cos(i)^2)*a_r^3*cos(2*aop + 2*f));
            
            % Continue with the rest
            n = sqrt(P.Con.primary.mu/a^3);
            [~,lpeSpec] = P.DynOeFourier([],IC,kMax);
            %             n = n + lpeSpec(11,1); % <-------------------  Work on this
            M = n*T+IC(6);
            k = 1:kMax;
            X = nan(6*P.Con.nSats,length(T));
            
            % Initial time
            trigMat = repmat([sin(k*M(1))./k/n;-cos(k*M(1))./k/n],6,1);
            trigsum1 = sum(lpeSpec(:,2:end).*trigMat,2);
            InitVal = [sum(trigsum1(1:2)); sum(trigsum1(3:4));...
                sum(trigsum1(5:6)); sum(trigsum1(7:8)); sum(trigsum1(9:10));...
                sum(trigsum1(11:12))];
            M2 = M;
            % Fix M
            
            for iTime = 1:length(T)
                trigMat = [sin(k*M(iTime))./k/n;-cos(k*M(iTime))./k/n];
                trigsum1 = sum(lpeSpec(11:12,2:end).*trigMat,2);
                trigsum2 = sum(trigsum1);
                M2(iTime) = lpeSpec(11,1)*T(iTime) + trigsum2 - InitVal(6) + M(iTime);
            end
            
            %
            for iTime = 1:length(T)
                trigMat = repmat([sin(k*M2(iTime))./k/n;-cos(k*M2(iTime))./k/n],6,1);
                trigsum1 = sum(lpeSpec(:,2:end).*trigMat,2);
                trigsum2 = [sum(trigsum1(1:2)); sum(trigsum1(3:4));...
                    sum(trigsum1(5:6)); sum(trigsum1(7:8)); sum(trigsum1(9:10));...
                    sum(trigsum1(11:12))];
                X(:,iTime) = IC + lpeSpec(1:2:11,1)*T(iTime) + trigsum2 - InitVal;
                %                 X(6,iTime) = X(6,iTime) + M2(iTime); % without M Fix
                X(6,iTime) = M2(iTime); % M fix
            end
            X(3:5,:) = wrapTo360(X(3:5,:)*180/pi);
            X(6,:) = X(6,:)*180/pi;
            X = X.';
            Time = T;
        end
        
        function [Time, X, hVec] = PropOeDeprit(P,nTime,nOrb)
            % Propagate a number of orbits with Deprit method.
            % Single Satellite only - could add more but probably not worth the time
            % ****** Has to start at perigee!!! ******
            mu = P.Con.primary.mu;
            Re = P.Con.primary.Re;
            J2 = P.Con.primary.J2;
            
            imagTol = 1e-12; % numerical tolerance for errors
            % Initial Conditions
            IC = P.Con.InitialOeOsc;
            sma = IC(1);
            ecc = IC(2);
            inc = IC(3);
            ran = IC(4);
            aop = IC(5);
%             man = IC(6);
%             f = me2ta(man,ecc);
            man = 0;
            f = 0;
            
            oeC = nan(nTime,6);
            oeW = nan(nTime,6);
            
            % Coordinate Switch
            radQ = sma*((1-ecc^2)/(1+ecc*cosd(f)));   % r
            aolQ = pi/180*(aop + f);                  % theta
            ranQ = pi/180*ran;                        % nu
            vraP = sqrt(mu/sma/(1-ecc^2))*ecc*sind(f);% R
            amoP = sqrt((1-ecc^2)*sma*mu);            % Theta
            amzP = amoP*cosd(inc);                    % N
            % Solution
            X = mu*J2*Re^2*(0.5-1.5*amzP^2/amoP^2);
            
            h = 0.5*vraP^2 + 0.5*amoP^2/radQ^2 - mu/radQ + ...
                0.25*mu*J2*Re^2/radQ^3*(1-3*amzP^2/amoP^2);
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
                if X < 0
                    k0 = (s2-s1)/(s3-s1);
                    n0 = (s2-s1)/s1;
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
                    n0 = (s3-s2)/s3;
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
                t = t.';
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
            if v0 ~= 0 % Check for singularity in unwrapping
                oeW(:,3) = v0/(2*pi)*unwrap(2*pi/v0*dRan.*signR)+ranQ;
            else % Polar orbit, RAAN is constant
                oeW(:,3) = ranQ;
            end
            oeW(:,4) = RVec;
            oeW(:,5) = amoP;
            oeW(:,6) = amzP;
            
            
            % Convert to Conventional
            oeC(:,1) = -mu*oeW(:,1).^2./(oeW(:,1).^2.*oeW(:,4).^2+...
                oeW(:,5).^2-2*mu*oeW(:,1));
            oeC(:,2) = sqrt(1-oeW(:,5).^2./(mu*oeC(:,1)));
            oeC(:,3) = acosd(oeW(:,6)./oeW(:,5));
            oeC(:,4) = wrapTo360(180/pi*oeW(:,3));
            oeC(:,5) = wrapTo360(180/pi*(oeW(:,2) - fVec));
            oeC(:,6) = 180/pi*unwrap(pi/180*ta2me(fVec*180/pi,oeC(:,2)));
            oeC(:,6) = oeC(:,6) - oeC(1,6); % force start at periapsis
            
            % Output
            Time = t;
            X = oeC;
            
        end
    end
    
    methods(Access = protected) % Equations of Motion
        % equations of motion functions should only be called by prop
        % functions, thus they are protected
        
        function dX = DynEciTb(P,t,X)  %#ok<INUSL>
            % move stuff around
            order = 6*P.Con.nSats;
            X2 = reshape(X,[6,P.Con.nSats]);
            Rv = [eye(3),zeros(3);
                zeros(3,6)]*X2;
            % get vector of R magnitudes
            r_vec = sqrt(dot(Rv,Rv,1));
            r_vec = reshape(repmat(r_vec,6,1),order,1);
            % move more stuff around
            R2 = repmat([1 1 1 0 0 0].',P.Con.nSats,1).*X;
            V2 = repmat([0 0 0 1 1 1].',P.Con.nSats,1).*X;
            % equation of motion
            dX = circshift(V2,-3) - circshift(P.Con.mu*R2./r_vec.^3,3);
        end
        
        function dX = DynEciJ2(P,t,X) %#ok<INUSL>
            % move stuff around
            order = 6*P.Con.nSats;
            X2 = reshape(X,[6,P.Con.nSats]);
            Rv = [eye(3),zeros(3);
                zeros(3,6)]*X2;
            % get vector of R magnitudes
            r = sqrt(dot(Rv,Rv,1));
            r = reshape(repmat(r,6,1),order,1);
            % move more stuff around
            R2 = repmat([1 1 1 0 0 0].',P.Con.nSats,1).*X;
            V2 = repmat([0 0 0 1 1 1].',P.Con.nSats,1).*X;
            Z  = repmat([0 0 1 0 0 0].',P.Con.nSats,1).*R2;
            Z2 = [zeros(3,2),ones(3,1),zeros(3);zeros(3,6)]*X2;
            Z2 = reshape(Z2,[order,1]);
            % equation of motion
            f_J2 = -P.Con.mu*P.Con.J2*P.Con.Re^2./r.^4.*...
                (3*Z./r + (-7.5*(Z2./r).^2 + 1.5).*R2./r);
            dX = circshift(V2,-3) + circshift(-P.Con.mu*R2./r.^3 + f_J2,3);
        end
        
        function dX = DynOeMean(P,t,X) %#ok<INUSL>
            % move stuff around
            order = 6*P.Con.nSats;
            X2 = reshape(X,[6,P.Con.nSats]);
            aV = X2(1,:);
            eV = X2(2,:);
            iV = X2(3,:);
            a = reshape(repmat(aV,6,1),order,1);
            e = reshape(repmat(eV,6,1),order,1);
            i = reshape(repmat(iV,6,1),order,1);
            % Derived values
            p = a.*(1-e.^2);
            n = sqrt(P.Con.mu./a.^3);
            eta = sqrt(1-e.^2);
            % Eq of motion
            dO = 180/pi*repmat([0,0,0,1,0,0].',P.Con.nSats,1).*...
                (-3/2*P.Con.J2.*(P.Con.Re./p).^2.*n.*cosd(i));
            dw = 180/pi*repmat([0,0,0,0,1,0].',P.Con.nSats,1).*...
                (3/4*P.Con.J2.*(P.Con.Re./p).^2.*n.*(5*cosd(i).^2-1));
            dM = 180/pi*repmat([0,0,0,0,0,1].',P.Con.nSats,1).*...
                (3/4*P.Con.J2.*(P.Con.Re./p).^2.*n.*eta.*(3*cosd(i).^2-1) + n);
            dX = dO + dw + dM;
        end
        
        function dX = DynOeOsc2(P,t,X)  %#ok<INUSL>
            OE  = reshape(X,6,P.Con.nSats);
            % Element vectors (angles already in radians)
            a   = OE(1,:);
            e   = OE(2,:);
            inc = OE(3,:);
            w   = OE(5,:);
            Me  = OE(6,:);
            th  = pi/180*me2ta(Me*180/pi,e);
            % Secondary definitions
            p = a.*(1-e.^2);
            h = sqrt(P.Con.mu*p);
            r = p./(1+e.*cos(th));
            MJ2 = -3*P.Con.mu./r.^4.*P.Con.J2.*P.Con.Re.^2;
            n = sqrt(P.Con.mu./a.^3);
            aol = th + w;
            % Forces
            f_r = MJ2/2.*(1-3*sin(inc).^2.*sin(aol).^2);
            f_th = MJ2.*sin(inc).^2.*sin(aol).*cos(aol);
            f_h = MJ2.*sin(inc).*cos(inc).*sin(aol);
            F = reshape([f_r;f_th;f_h],3*P.Con.nSats,1);
            % GVE Sparsified Matrix
            B11 = vec2sparse(2*a.^2./h.*e.*sin(th),[6,3],[1,1]);
            B12 = vec2sparse(2*a.^2./h.*(1 + e.*cos(th)),[6,3],[1,2]);
            B21 = vec2sparse(p./h.*sin(th),[6,3],[2,1]);
            B22 = vec2sparse(r./h.*(e + 2*cos(th) + e.*cos(th).^2),[6,3],[2,2]);
            B33 = vec2sparse(r./h.*cos(aol),[6,3],[3,3]);
            B43 = vec2sparse(r.*sin(aol)./(h.*sin(inc)),[6,3],[4,3]);
            B51 = vec2sparse(-p./(h.*e).*cos(th),[6,3],[5,1]);
            B52 = vec2sparse(r./(h.*e).*(2+e.*cos(th)).*sin(th),[6,3],[5,2]);
            B53 = vec2sparse(-r./h.*sin(aol).*cos(inc)./sin(inc),[6,3],[5,3]);
            B61 = vec2sparse((p.*cos(th)-2*r.*e)./(n.*a.^2.*e),[6,3],[6,1]);
            B62 = vec2sparse(-(p+r).*sin(th)./(n.*a.^2.*e),[6,3],[6,2]);
            
            B = B11 + B12 + B21 + B22 + B33 + B43 + B51 + B52 + B53 + B61 + B62;
            
            % Kepler Solution
            K = reshape([zeros(5,P.Con.nSats);n],6*P.Con.nSats,1);
            
            % Equations of Motion
            dX = B*F + K;
        end
        
        function dX = DynOeOsc(P,t,X)  %#ok<INUSL>
            OE  = reshape(X,6,P.Con.nSats);
            % Element vectors (angles already in radians)
            a   = OE(1,:);
            e   = OE(2,:);
            inc = OE(3,:);
            w   = OE(5,:);
            Me  = OE(6,:);
            th  = pi/180*me2ta(Me*180/pi,e);
            % Secondary definitions
            p = a.*(1-e.^2);
            h = sqrt(P.Con.mu*p);
            r = p./(1+e.*cos(th));
            MJ2 = -3*P.Con.mu./r.^4.*P.Con.J2.*P.Con.Re.^2;
            n = sqrt(P.Con.mu./a.^3);
            aol = th + w;
            % Forces
            fR = MJ2/2.*(1-3*sin(inc).^2.*sin(aol).^2);
            fTh = MJ2.*sin(inc).^2.*sin(aol).*cos(aol);
            fH = MJ2.*sin(inc).*cos(inc).*sin(aol);
            % Element Rates
            da = 2*a.^2./h.*e.*sin(th).*fR +...
                2*a.^2./h.*(1 + e.*cos(th)).*fTh;
            de = p./h.*sin(th).*fR +...
                r./h.*(e + 2*cos(th) + e.*cos(th).^2).*fTh;
            di = r./h.*cos(aol).*fH;
            dO = r.*sin(aol)./(h.*sin(inc)).*fH;
            dw = -p./(h.*e).*cos(th).*fR +...
                r./(h.*e).*(2+e.*cos(th)).*sin(th).*fTh +...
                -r./h.*sin(aol).*cos(inc)./sin(inc).*fH;
            dM = (p.*cos(th)-2*r.*e)./(n.*a.^2.*e).*fR +...
                -(p+r).*sin(th)./(n.*a.^2.*e).*fTh;
            
            dOe = reshape([da;de;di;dO;dw;dM + n],6*P.Con.nSats,1);
            % Equations of Motion
            dX = dOe;
        end
        
        function dX = DynEciJ3(P,t,X) %#ok<INUSL>
            % move stuff around
            order = 6*P.Con.nSats;
            X2 = reshape(X,[6,P.Con.nSats]);
            Rv = [eye(3),zeros(3);
                zeros(3,6)]*X2;
            % get vector of R magnitudes
            r = sqrt(dot(Rv,Rv,1));
            r = reshape(repmat(r,6,1),order,1);
            % move more stuff around
            R2 = repmat([1 1 1 0 0 0].',P.Con.nSats,1).*X;
            V2 = repmat([0 0 0 1 1 1].',P.Con.nSats,1).*X;
            Z  = repmat([0 0 1 0 0 0].',P.Con.nSats,1).*R2;
            Z2 = [zeros(3,2),ones(3,1),zeros(3);zeros(3,6)]*X2;
            Z2 = reshape(Z2,[order,1]);
            % equation of motion
            f_J2 = -P.Con.mu*P.Con.J2*P.Con.Re^2./r.^4.*...
                (3*Z./r + (-7.5*(Z2./r).^2 + 1.5).*R2./r);
            
            f_J3 = -5/2*P.Con.mu*P.Con.J3*P.Con.Re^3./r.^7.*...
                ((3*Z2 - 7*Z2.^2./r.^2).*R2 + ...
                (3 - 3/5*r.^2./(Z2.^2)).*Z.^2);
            
            dX = circshift(V2,-3) + ...
                circshift(-P.Con.mu*R2./r.^3 + f_J2 + f_J3,3);
        end
        
        function [dX,lpeSpec] = DynOeFourier(P,t,X,kMax) %#ok<INUSL>
            
            % handle elements vector
            a = X(1);
            e = X(2);
            i = X(3);
            aop = X(5);
            M = X(6);
            b = (1-sqrt(1-e^2))/e;
            
            % constant potential values
            R = -P.Con.primary.mu*P.Con.primary.J2*P.Con.primary.Re^2/2/a^3; % common factor
            dRda = -3*R/a;
            
            R0 = -(3*cos(i)^2-1)/2/(1-e^2)^(3/2); % 0 freq element
            dR0di = 3*cos(i)*sin(i)/(1-e^2)^(3/2);
            dR0de = -3*e*(3*cos(i)^2-1)/2/(1-e^2)^(5/2); %
            dR0do = 0;
            dR0dl = 0;
            
            % constant vectors
            m2 = (0:4).';
            m3 = (0:6).';
            m4 = (0:8).';
            m5 = (0:10).';
            
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
            
            dC1xA1di = 36*sin(i)*cos(i)*cos(aop)^2/(1-e^2)^2;
            dC1xA1do = -36*sin(i)^2*sin(aop)*cos(aop)/(1-e^2)^2;
            dC1xA1de =  24*e*(3*sin(i)^2*cos(aop)^2-1)/(1-e^2)^3;
            
            C = [(9*e^2*sin(i)^2*cos(aop)^2+3*sin(i)^2*(sin(aop)^2-cos(aop)^2)-3*e^2)/2/(1-e^2)^(7/2);
                (3*e^3*sin(i)^2*cos(aop)^2+9*e*sin(i)^2*(sin(aop)^2-cos(aop)^2)-e^3)/4/(1-e^2)^4;
                (9*e^2*sin(i)^2*(sin(aop)^2-cos(aop)^2))/8/(1-e^2)^(9/2);
                (3*e^3*sin(i)^2*(sin(aop)^2-cos(aop)^2))/16/(1-e^2)^5];
            
            dCdi = [3*(3*e^2*sin(i)*cos(i)*cos(aop)^2+sin(i)*cos(i)*(sin(aop)^2-cos(aop)^2))/(1-e^2)^(7/2);
                3*(e^3*sin(i)*cos(i)*cos(aop)^2+3*e*sin(i)*cos(i)*(sin(aop)^2-cos(aop)^2))/2/(1-e^2)^4;
                9*e^2*sin(i)*cos(i)*(sin(aop)^2-cos(aop)^2)/4/(1-e^2)^(9/2);
                3*e^3*sin(i)*cos(i)*(sin(aop)^2-cos(aop)^2)/8/(1-e^2)^5];
            
            dCdo = [3*(-3*e^2+2)*sin(i)^2*sin(aop)*cos(aop)/(1-e^2)^(7/2);
                -3*e*(e^2-6)*sin(i)^2*sin(aop)*cos(aop)/2/(1-e^2)^4;
                9*e^2*sin(i)^2*sin(aop)*cos(aop)/2/(1-e^2)^(9/2);
                3*e^3*sin(i)^2*sin(aop)*cos(aop)/4/(1-e^2)^5];
            
            dCde = [3*e*(15*e^2*sin(i)^2*cos(aop)^2 +7*sin(i)^2*sin(aop)^2 -sin(i)^2*cos(aop)^2 -5*e^2 -2)/2/(1-e^2)^(9/2);
                ((15*e^4-54*e^2-9)*sin(i)^2*cos(aop)^2 + (63*e^2+9)*sin(i)^2*sin(aop)^2 -5*e^4 -3*e^2)/4/(1-e^2)^5;
                -9*e*sin(i)^2*(cos(aop)^2-sin(aop)^2)*(7*e^2+2)/8/(1-e^2)^(11/2);
                -3*e^2*sin(i)^2*(cos(aop)^2-sin(aop)^2)*(7*e^2+3)/16/(1-e^2)^6];
            
            
            S = -[6*sin(i)^2*sin(aop)*cos(aop)/2/(1-e^2)^3;
                18*e*sin(i)^2*sin(aop)*cos(aop)/4/(1-e^2)^(7/2);
                18*e^2*sin(i)^2*sin(aop)*cos(aop)/8/(1-e^2)^4;
                6*e^3*sin(i)^2*sin(aop)*cos(aop)/16/(1-e^2)^(9/2)];
            
            dSdi = -[12*sin(i)*cos(i)*sin(aop)*cos(aop)/2/(1-e^2)^3;
                36*e*sin(i)*cos(i)*sin(aop)*cos(aop)/4/(1-e^2)^(7/2);
                36*e^2*sin(i)*cos(i)*sin(aop)*cos(aop)/8/(1-e^2)^4;
                12*e^3*sin(i)*cos(i)*sin(aop)*cos(aop)/16/(1-e^2)^(9/2)];
            
            dSdo = -[6*sin(i)^2*(2*cos(aop)^2-1)/2/(1-e^2)^3;
                18*e*sin(i)^2*(2*cos(aop)^2-1)/4/(1-e^2)^(7/2);
                18*e^2*sin(i)^2*(2*cos(aop)^2-1)/8/(1-e^2)^4;
                6*e^3*sin(i)^2*(2*cos(aop)^2-1)/16/(1-e^2)^(9/2)];
            
            dSde = [-18*e*sin(i)^2*sin(aop)*cos(aop)/(1-e^2)^4;
                -9*sin(i)^2*sin(aop)*cos(aop)*(6*e^2+1)/2/(1-e^2)^(9/2);
                -9*e*sin(i)^2*sin(aop)*cos(aop)*(3*e^2+1)/2/(1-e^2)^5;
                -9*e^2*sin(i)^2*sin(aop)*cos(aop)*(2*e^2+1)/8/(1-e^2)^(11/2)];
            
            dJ2daFreq = dRda*[[R0;0],zeros(2,kMax)];
            dJ2deFreq = R*[[dR0de;0],zeros(2,kMax)];
            dJ2diFreq = R*[[dR0di;0],zeros(2,kMax)];
            dJ2doFreq = R*[[dR0do;0],zeros(2,kMax)];
            dJ2dlFreq = R*[[dR0dl;0],zeros(2,kMax)];
            
            k = 1;
            
            while k <= kMax
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
                
                Jk = besselj(k,k*e);
                dJkde = 0.5*(besselj(k-1,k*e) - besselj(k+1,k*e));
                Jn = besselj(n,-k*e);
                dJnde = 0.5*(-besselj(1-n,-k*e) - besselj(n+1,-k*e)); % uses identity only correct for n = 0
                
                Akda = C1xA1*Jk +...
                    Jn*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
                Akde = (C1xA1*k*dJkde + dC1xA1de*Jk) + ...
                    -k*dJnde*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]+...
                    Jn*(dCde.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] + ...
                    C.'*[da2de.'*g2+a2.'*dg2de; da3de.'*g3+a3.'*dg3de; da4de.'*g4+a4.'*dg4de; da5de.'*g5+a5.'*dg5de]);
                Akdi = dC1xA1di*Jk +...
                    Jn*dCdi.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
                Akdo = dC1xA1do*Jk +...
                    Jn*dCdo.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
                Akdl = Akda;
                
                Bkda = Jn*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                Bkde = -k*dJnde*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5]+...
                    Jn*(dSde.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] + ...
                    S.'*[db1de.'*g2+b1.'*dg2de; db2de.'*g3+b2.'*dg3de; db3de.'*g4+b3.'*dg4de; db4de.'*g5+b4.'*dg5de]);
                Bkdi = Jn*dSdi.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                Bkdo = Jn*dSdo.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                Bkdl = Bkda;
                
                n = 1;
                while n <= k + 8
                    % positive n
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
                    
                    Jn = besselj(n,-k*e);
                    dJnde = 0.5*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
                    
                    dAkda = Jn*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
                    dAkde = -k*dJnde*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]+...
                        Jn*(dCde.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] + ...
                        C.'*[da2de.'*g2+a2.'*dg2de; da3de.'*g3+a3.'*dg3de; da4de.'*g4+a4.'*dg4de; da5de.'*g5+a5.'*dg5de]);
                    dAkdi = Jn*dCdi.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
                    dAkdo = Jn*dCdo.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
                    
                    
                    dBkda = Jn*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBkde = -k*dJnde*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5]+...
                        Jn*(dSde.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] + ...
                        S.'*[db1de.'*g2+b1.'*dg2de; db2de.'*g3+b2.'*dg3de; db3de.'*g4+b3.'*dg4de; db4de.'*g5+b4.'*dg5de]);
                    dBkdi = Jn*dSdi.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBkdo = Jn*dSdo.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    
                    
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
                    
                    dAkda = dAkda + Jn*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
                    dAkde = dAkde + k*dJnde*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]+...
                        Jn*(dCde.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] + ...
                        C.'*[da2de.'*g2+a2.'*dg2de; da3de.'*g3+a3.'*dg3de; da4de.'*g4+a4.'*dg4de; da5de.'*g5+a5.'*dg5de]);
                    dAkdi = dAkdi + Jn*dCdi.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
                    dAkdo = dAkdo + Jn*dCdo.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
                    dAkdl = dAkda;
                    
                    dBkda = dBkda + Jn*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBkde = dBkde + k*dJnde*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5]+...
                        Jn*(dSde.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] + ...
                        S.'*[db1de.'*g2+b1.'*dg2de; db2de.'*g3+b2.'*dg3de; db3de.'*g4+b3.'*dg4de; db4de.'*g5+b4.'*dg5de]);
                    dBkdi = dBkdi + Jn*dSdi.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBkdo = dBkdo + Jn*dSdo.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBkdl = dBkda;
                    
                    Akda = Akda + dAkda;
                    Akde = Akde + dAkde;
                    Akdi = Akdi + dAkdi;
                    Akdo = Akdo + dAkdo;
                    Akdl = Akdl + dAkdl;
                    
                    Bkda = Bkda + dBkda;
                    Bkde = Bkde + dBkde;
                    Bkdi = Bkdi + dBkdi;
                    Bkdo = Bkdo + dBkdo;
                    Bkdl = Bkdl + dBkdl;
                    
                    n = n+1;
                end
                dJ2daFreq(:,k+1) = dRda*[Akda;Bkda];
                dJ2deFreq(:,k+1) = R*[Akde;Bkde];
                dJ2diFreq(:,k+1) = R*[Akdi;Bkdi];
                dJ2doFreq(:,k+1) = R*[Akdo;Bkdo];
                dJ2dlFreq(:,k+1) = R*[k*Bkdl;-k*Akdl];
                
                k = k+1;
            end
            k = 0:kMax;
            trigMat = [cos(k*M);sin(k*M)];
            
            dJ2da = sum(dJ2daFreq.*trigMat,'all');
            dJ2de = sum(dJ2deFreq.*trigMat,'all');
            dJ2di = sum(dJ2diFreq.*trigMat,'all');
            dJ2do = sum(dJ2doFreq.*trigMat,'all');
            dJ2dl = sum(dJ2dlFreq.*trigMat,'all');
            
            dX = zeros(6,1);
            n = sqrt(P.Con.primary.mu/a^3);
            eta = sqrt(1-e^2);
            
            dX(1) = 2/n/a*dJ2dl;
            dX(2) = eta^2/n/a^2/e*dJ2dl - eta/n/a^2/e*dJ2do;
            dX(3) = cos(i)/n/a^2/eta/sin(i)*dJ2do;
            dX(4) = 1/n/a^2/eta/sin(i)*dJ2di;
            dX(5) = eta/n/a^2/e*dJ2de - cos(i)/n/a^2/eta/sin(i)*dJ2di;
            dX(6) = n -2/n/a*dJ2da -eta^2/n/a^2/e*dJ2de;
            
            % Fix a,n - Not clear if helpful
            %             f = me2ta(M,e);
            %             a_r = (1+e.*cos(f))./eta.^2;
            %             g2 = -P.Con.primary.J2/2*(P.Con.primary.Re/a)^2;
            %             a = a + a*g2*((3*cos(i)^2-1).*(a_r^3 - 1/eta^3) ...
            %                 + 3*(1-cos(i)^2)*a_r^3*cos(2*aop + 2*f));
            %             n = sqrt(P.Con.primary.mu/a^3);
            
            lpeSpec = nan(12,kMax+1);
            
            lpeSpec(1:2,:) = 2/n/a*dJ2dlFreq;
            lpeSpec(3:4,:) = eta^2/n/a^2/e*dJ2dlFreq - eta/n/a^2/e*dJ2doFreq;
            lpeSpec(5:6,:) = cos(i)/n/a^2/eta/sin(i)*dJ2doFreq;
            lpeSpec(7:8,:) = 1/n/a^2/eta/sin(i)*dJ2diFreq;
            lpeSpec(9:10,:) = eta/n/a^2/e*dJ2deFreq - cos(i)/n/a^2/eta/sin(i)*dJ2diFreq;
            lpeSpec(11:12,:) = -2/n/a*dJ2daFreq -eta^2/n/a^2/e*dJ2deFreq;
        end
        
    end
    
    methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(P)  %#ok<MANU>
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Propagation Definitions';
            propgroups(1).PropertyList = {'relTol','absTol','Con'};
        end
    end
end