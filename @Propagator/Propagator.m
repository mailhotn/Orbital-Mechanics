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
                    relTol     = 1e-11;
                    absTol     = 1e-12;
                    Constellation = WalkerConstellation;
                case 1 % default tolerances
                    relTol     = 1e-11;
                    absTol     = 1e-12;
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
            % Normalize
            dScale = P.Con.Re;
            tScale = sqrt(P.Con.Re^3/P.Con.mu);
            eciScale = repmat([dScale*ones(3,1);dScale/tScale*ones(3,1)],P.Con.nSats,1);

            IC = reshape(P.Con.InitialStateEci,[6*P.Con.nSats,1])./eciScale;
            T = T/tScale;
            % Prop Normalized
            [Time, X] = ode78(@P.DynEciJ2,T,IC,opts);
            % De-Normalize
            Time = Time*tScale;
            X = X.*eciScale.';
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
        
        function [Time, X] = PropOeMeanShort(P,T)
            % Directly calculates linear solution without numerical
            % integration, using Kozai w/o Long Period. Avoids singularity
            IC = reshape(P.Con.InitialOeMeanShort,[6*P.Con.nSats,1]);
            % move stuff around
            order = 6*P.Con.nSats;
            X2 = reshape(IC,[6,P.Con.nSats]);
            aV = X2(1,:);
            eV = X2(2,:);
            iV = X2(3,:);
            sma = reshape(repmat(aV,6,1),order,1);
            ecc = reshape(repmat(eV,6,1),order,1);
            inc = reshape(repmat(iV,6,1),order,1);
            % Derived values
            p = sma.*(1-ecc.^2);
            n = sqrt(P.Con.mu./sma.^3);
            eta = sqrt(1-ecc.^2);
            % Eq of motion
            dO = 180/pi*repmat([0,0,0,1,0,0].',P.Con.nSats,1).*...
                (-3/2*P.Con.J2.*(P.Con.Re./p).^2.*n.*cosd(inc));
            dw = 180/pi*repmat([0,0,0,0,1,0].',P.Con.nSats,1).*...
                (3/4*P.Con.J2.*(P.Con.Re./p).^2.*n.*(5*cosd(inc).^2-1));
            dM = 180/pi*repmat([0,0,0,0,0,1].',P.Con.nSats,1).*...
                (3/4*P.Con.J2.*(P.Con.Re./p).^2.*n.*eta.*(3*cosd(inc).^2-1) + n);
            dX = dO + dw + dM;
            % "Propagate"
            X = zeros(order,length(T));
            X(:,1) = IC;
            X(:,2:end) = X(:,1) + dX*T(2:end);
            X = X.';
            Time = T;
        end
        
        function [Time, X] = PropOeOsc(P,T)
            % Numerically propagate GVE for osculating elements
            % Singular in e, not i
            % Use this I think
            opts = odeset('reltol',P.relTol,'abstol',P.absTol);
            % Normalize
            dScale = P.Con.Re;
            tScale = sqrt(P.Con.Re^3/P.Con.mu);
            oeScale = repmat([dScale;ones(5,1)],P.Con.nSats,1);

            OE = P.Con.InitialOeOsc;
            OE(3:end,:) = OE(3:end,:)*pi/180;
            IC = reshape(OE,[6*P.Con.nSats,1])./oeScale;
            T = T/tScale;
            % Prop Normalized
            [Time, X] = ode78(@P.DynOeOsc,T,IC,opts);
            inddeg = reshape((3:6).'+(0:P.Con.nSats-1)*6,4*P.Con.nSats,1);
            indWrap = reshape((3:5).'+(0:P.Con.nSats-1)*6,3*P.Con.nSats,1);
            X(:,inddeg) = (180/pi*X(:,inddeg));
            X(:,indWrap) = wrapTo360(X(:,indWrap));
            % De-Normalize
            Time = Time*tScale;
            X = X.*oeScale.';
        end
        
        function [Time, X] = PropOeOsc2(P,T)
            % Numerically propagate GVE for osculating elements
            % Tried different way of writing equations, didnt work well
            % probably worse?
            opts = odeset('reltol',P.relTol,'abstol',P.absTol);
            OE = P.Con.InitialOeOsc;
            OE(3:end,:) = OE(3:end,:)*pi/180;
            IC = reshape(OE,[6*P.Con.nSats,1]);
            [Time, X] = ode78(@P.DynOeOsc2,T,IC,opts);
            inddeg = reshape((3:6).'+(0:P.Con.nSats-1)*6,4*P.Con.nSats,1);
            X(:,inddeg) = wrapTo360(180/pi*X(:,inddeg));
        end
        
        function [Time, X] = PropOeOsc3(P,T)
            % Numerically propagate LPE for osculating elements
            % Singular in e, not i
            % Use this I think
            opts = odeset('reltol',P.relTol,'abstol',P.absTol);
            % Normalize
            dScale = P.Con.Re;
            tScale = sqrt(P.Con.Re^3/P.Con.mu);
            oeScale = repmat([dScale;ones(5,1)],P.Con.nSats,1);

            OE = P.Con.InitialOeOsc;
            OE(3:end,:) = OE(3:end,:)*pi/180;
            IC = reshape(OE,[6*P.Con.nSats,1])./oeScale;
            T = T/tScale;
            % Prop Normalized
            [Time, X] = ode78(@P.DynOeOsc3,T,IC,opts);
            inddeg = reshape((3:6).'+(0:P.Con.nSats-1)*6,4*P.Con.nSats,1);
            indWrap = reshape((3:5).'+(0:P.Con.nSats-1)*6,3*P.Con.nSats,1);
            X(:,inddeg) = (180/pi*X(:,inddeg));
            X(:,indWrap) = wrapTo360(X(:,indWrap));
            % De-Normalize
            Time = Time*tScale;
            X = X.*oeScale.';
        end
        
        function [Time, X, PNS] = PropOePns(P,T)
            % Numerically propagate GVE for oscullating elements
            % Uses non-singular Polar-Nodal elements
            opts = odeset('reltol',P.relTol,'abstol',P.absTol);
            % Normalize
            dScale = P.Con.Re;
            tScale = sqrt(P.Con.Re^3/P.Con.mu);
            pnsScale = repmat([dScale;1;dScale/sqrt(tScale);...
                dScale/tScale;dScale^2/tScale;dScale/sqrt(tScale)],P.Con.nSats,1);

            OE = P.Con.InitialOeOsc;
            OE(6,:) = me2ta(OE(6,:),OE(2,:));
            PNS = oe2pns(OE,P.Con.primary);
            IC = reshape(PNS,[6*P.Con.nSats,1])./pnsScale;
            T = T/tScale;
            [Time, X] = ode78(@P.DynOePns,T,IC,opts);
            Time = Time*tScale;
            X = X.*pnsScale.';
            PNS = X;
            X = reshape(X.',6, P.Con.nSats*length(Time));
            X = pns2oe(X,P.Con.primary);
            X = ta2me(X);
            X = reshape(X,6*P.Con.nSats,length(Time)).';
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
            %             X = X;
            %             X(3:end,:) = wrapTo360(X(3:end,:));
        end
        
        function [Time, X] = PropOeFourier2(P,T,kMax)
            % Refactoring Bullshit
            [Time, X] = PropOeFourier(P,T,kMax);
        end
        
        function [Time, X] = PropOeFourier(P,T,kMax)
            % Propagate for time T using Fourier LPE
            % Assume constant coefficients
            
            % Handle Initial conditions
            ic = reshape(P.Con.InitialOeOsc,[6*P.Con.nSats,1]);
            icOsc = ic;
            icOsc(3:end) = icOsc(3:end)*pi/180;
            ic = osc2meSP(ic);
            ic(3:end) = ic(3:end)*pi/180;
            
            a = ic(1);
            
            %             % Get mean a by substracting short period variations. The mean
            %             % mean motion resulting *should* be the correct frequency for M
            %             e = ic(2);
            %             i = ic(3);
            %             aop = ic(5);
            %
            %             eta = sqrt(1-e^2);
            %             f = pi/180*me2ta(ic(6)*180/pi,e);
            %             a_r = (1+e.*cos(f))./eta.^2;
            %             g2 = -P.Con.primary.J2/2*(P.Con.primary.Re/a)^2;
            % %             a = a + 3*g2*a*(1-1.5*sin(i)^2)/eta^3; % Kozai - weird
            %             a = a + a*g2*((3*cos(i)^2-1).*(a_r^3 - 1/eta^3) ...
            %                 + 3*(1-cos(i)^2)*a_r^3*cos(2*aop + 2*f));
            %             ic(1) = a; % Reassign a given averaging
            
            % Continue with the rest
            n = sqrt(P.Con.primary.mu/a^3);
            %             n = sqrt(P.Con.primary.mu/a^3)*(1+3*g2*(1-1.5*sin(i)^2)/eta^3); % kozai Fix
            [~,lpeSpec] = P.DynOeFourier([],ic,kMax);
            %                         n = n + lpeSpec(11,1); % <-------------------  Work on this
            M = n*T+icOsc(6);
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
            Sk = sin(k.'*M)./k.'/n;
            Ck = -cos(k.'*M)./k.'/n;
            AkM = lpeSpec(11,2:end);
            BkM = lpeSpec(12,2:end);
            M2 = lpeSpec(11,1)*T + AkM*Sk + BkM*Ck - InitVal(6) + M;
            
            % Calculate all elements
            Sk = sin(k.'*M2)./k.'/n;
            Ck = -cos(k.'*M2)./k.'/n;
            Ak = lpeSpec(1:2:11,2:end);
            Bk = lpeSpec(2:2:12,2:end);
            
            X = icOsc + lpeSpec(1:2:11,1)*T + Ak*Sk + Bk*Ck -InitVal;
            X(6,:) = M2;
            
            
            X(3:5,:) = wrapTo360(X(3:5,:)*180/pi);
            X(6,:) = X(6,:)*180/pi;
            X = X.';
            Time = T;
        end
        
        function [Time, X] = PropOeFourier2Ord(P,T,kMax)
            % Propagate for time T using Fourier LPE
            % Assume coefficients varying in AOP up to first order

            %% Handle Initial conditions
            icVec = reshape(P.Con.InitialOeOsc,[6*P.Con.nSats,1]); % ic of all sats
            X = nan(6*P.Con.nSats,length(T));
            % for satellites - inefficient, but no good option
            for iSat = 1:P.Con.nSats
                icOsc = icVec((1:6)+6*(iSat-1));

                [freq0,lpeSpec] = P.DynOeFourier2Ord(T,icOsc,kMax);

                % Average out SP from initial conditions
                icM = osc2meSP(icOsc); % change to numerical mean
                icOsc(3:end) = icOsc(3:end)*pi/180;
                icM(3:end) = icM(3:end)*pi/180;

                smaM = icM(1);

                % Continue with the rest
                nM = sqrt(P.Con.primary.mu/smaM^3);

                M = nM*T+icOsc(6);
                k = 1:kMax;
                Xi = nan(6,length(T));

                % Initial M
                trigMat = [sin(k.'*M(1))./k.'/nM;-cos(k.'*M(1))./k.'/nM];
                InitM = sum(lpeSpec((10*kMax+1):end,1).*trigMat);


                % Fix M
                Sk = sin(k.'*M)./k.'/nM;
                Ck = -cos(k.'*M)./k.'/nM;
                AkM = lpeSpec((10*kMax+1):11*kMax,:);
                BkM = lpeSpec((11*kMax+1):12*kMax,:);
                M2 = freq0(6)*T + sum(AkM.*Sk + BkM.*Ck) - InitM + M;

                % Initial time - With fixed M
                trigMat = repmat([sin(k.'*M2(1))./k.'/nM;-cos(k.'*M2(1))./k.'/nM],6,1);
                trigsum1 = lpeSpec(:,1).*trigMat;
                % CHANGE TO MATRIX
                % MULTIPLICATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                InitVal = [sum(trigsum1(1:2*kMax)); sum(trigsum1((2*kMax+1):4*kMax));...
                    sum(trigsum1((4*kMax+1):6*kMax)); sum(trigsum1((6*kMax+1):8*kMax));...
                    sum(trigsum1((8*kMax+1):10*kMax)); sum(trigsum1((10*kMax+1):12*kMax))];


                % Calculate all elements
                Sk = repmat(sin(k.'*M2)./k.'/nM,6,1);
                Ck = repmat(-cos(k.'*M2)./k.'/nM,6,1);
                iAk = (1:(6*kMax)) + reshape(repmat((0:5)*kMax,kMax,1),1,kMax*6); % Indices for Ak in Spectrum
                iBk = (1:(6*kMax)) + reshape(repmat((1:6)*kMax,kMax,1),1,kMax*6);
                Ak = lpeSpec(iAk,:);
                Bk = lpeSpec(iBk,:);
                sumMat = blkdiag(ones(1,kMax),ones(1,kMax),ones(1,kMax),ones(1,kMax),ones(1,kMax),ones(1,kMax));

                Xi = icOsc + freq0*T + sumMat*(Ak.*Sk + Bk.*Ck) - InitVal;
                %             X(6,:) = M2;
                Xi(6,:) = Xi(6,:) + M - icOsc(6); % subtract M(0)? test with nonzero IC


                Xi(3:5,:) = wrapTo360(Xi(3:5,:)*180/pi);
                Xi(6,:) = Xi(6,:)*180/pi;
                % finished Xi, assign then iterate
                X((1:6)+6*(iSat-1),:) = Xi; 
            end
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
                    % Precalc Elliptics
                    [eK,eE,eP] = elliptic123(k0,-n0); % complete
                    [eFi,eEi,ePi] = elliptic123(phi,k0,-n0); % incomplete
                    %                     eFi = ellipticF(phi,k0);
                    %                     eEi = ellipticE(phi,k0);
                    
                    % RAAN solution
                    Iv0 = 2*sqrt(s3/-X)*(sqrt(s3/(s3-s1))*eK-...
                        sqrt((s3-s1)/s3)*eE);
                    Iv = Iv0 - 2*sqrt(s3/-X)*(sqrt(s3/(s3-s1))*eFi-...
                        sqrt((s3-s1)/s3)*eEi);
                    v0 = -3/2*mu*J2*Re^2*amzP/amoP^2*Iv0;
                    
                    % AOL solution
                    Ith0 = 2/(sqrt(-X)*sqrt(s3-s1))*eK;
                    Ith = Ith0 - 2/(sqrt(-X)*sqrt(s3-s1))*eFi;
                    th0 = amoP*Ith0 - v0*amzP/amoP;
                    
                    % Time Solution
                    T0 = 1/(sqrt(-X)*s1^2*sqrt(s2-s1))*sqrt(k0)/(1+n0)*...
                        (((3+2*n0)*k0+(2+n0)*n0)/(k0+n0)*eP + ...
                        n0/(k0+n0)*eE - eK);
                    It = T0 - 1/(sqrt(-X)*s1^2*sqrt(s2-s1))*sqrt(k0)/(1+n0)*...
                        (n0/2*n0/(k0+n0)*sqrt(1-k0*sin(phi).^2)./(1+n0*sin(phi).^2).*sin(2*phi)...
                        +((3+2*n0)*k0+(2+n0)*n0)/(k0+n0)*ePi ...
                        + n0/(k0+n0)*eEi -eFi);
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
                    % Precalc Elliptics
                    [eK,eE,eP] = elliptic123(k0,n0);
                    [eFi,eEi,ePi] = elliptic123(phi,k0,n0);
                    %                     eFi = ellipticF(phi,k0);
                    %                     eEi = ellipticE(phi,k0);
                    
                    % RAAN solution
                    Iv0 = 2*sqrt(s1/X)*(sqrt(s1/(s3-s1))*eK-...
                        sqrt((s3-s1)/s1)*eE);
                    Iv = 2*sqrt(s1/X)*(sqrt(s1/(s3-s1))*eFi-...
                        sqrt((s3-s1)/s1)*eEi);
                    v0 = -3/2*mu*J2*Re^2*amzP/amoP^2*Iv0; %v*
                    
                    % AOL solution
                    Ith0 = 2/(sqrt(X)*sqrt(s3-s1))*eK;
                    Ith = 2/(sqrt(X)*sqrt(s3-s1))*eFi;
                    th0 = amoP*Ith0 - v0*amzP/amoP; %th*
                    
                    % Time Solution
                    T0 = 1/(sqrt(X)*s3^2*sqrt(s3-s2))*sqrt(k0)/(1-n0)*...
                        (((3-2*n0)*k0-(2-n0)*n0)/(k0-n0)*eP - ...
                        n0/(k0-n0)*eE - eK);
                    It = 1/(sqrt(X)*s3^2*sqrt(s3-s2))*sqrt(k0)/(1-n0)*...
                        (n0/2*n0/(k0-n0)*sqrt(1-k0*sin(phi).^2)./(1-n0*sin(phi).^2).*sin(2*phi)...
                        +(3*k0-2*n0-(2*k0-n0)*n0)/(k0-n0)*ePi ...
                        - n0/(k0-n0)*eEi -eFi);
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
            
            
            % Convert to Conventional elements
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
            % Normalized Parameters
            mu = 1;
            Re = 1;
            % move stuff around
            order = 6*P.Con.nSats;
            X2 = reshape(X,[6,P.Con.nSats]);
            Rv = [eye(3),zeros(3);
                zeros(3,6)]*X2;
            % get vector of R magnitudes
%             r = sqrt(dot(Rv,Rv,1));
            r = vecnorm(Rv);
            r = reshape(repmat(r,6,1),order,1);
            % move more stuff around
            R2 = reshape(Rv,order,1);
            V2 = reshape([zeros(3,6);zeros(3),eye(3)]*X2,order,1);
            Z  = repmat([0 0 1 0 0 0].',P.Con.nSats,1).*R2;
            Z2 = reshape([zeros(3,2),ones(3,1),zeros(3);zeros(3,6)]*X2,order,1);
            % equation of motion
            f_J2 = -mu*P.Con.J2*Re^2./r.^4.*...
                (3*Z./r + (-7.5*(Z2./r).^2 + 1.5).*R2./r);
            dX = circshift(V2,-3) + circshift(-mu*R2./r.^3 + f_J2,3);
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
        
        function dX = DynOeOsc(P,t,X)  %#ok<INUSL>
            OE  = reshape(X,6,P.Con.nSats);
            % Normalized Parameters
            mu = 1;
            Re = 1;
            % Element vectors (angles already in radians)
            a   = OE(1,:);
            e   = OE(2,:);
            inc = OE(3,:);
            w   = OE(5,:);
            Me  = OE(6,:);
            th  = pi/180*me2ta(Me*180/pi,e);
            % Secondary definitions
            p = a.*(1-e.^2);
            h = sqrt(mu*p);
            r = p./(1+e.*cos(th));
            MJ2 = -3*mu./r.^4.*P.Con.J2.*Re.^2;
            n = sqrt(mu./a.^3);
            aol = th + w;
            % Forces
            fR = MJ2/2.*(1-3*sin(inc).^2.*sin(aol).^2);
            fTh = MJ2.*sin(inc).^2.*sin(aol).*cos(aol);
            fH_si = MJ2.*cos(inc).*sin(aol); % sin(inc) taken out to avoid singularity
            % Element Rates
            da = 2*a.^2./h.*e.*sin(th).*fR +...
                2*a.^2./h.*(1 + e.*cos(th)).*fTh;
            de = p./h.*sin(th).*fR +...
                r./h.*(e + 2*cos(th) + e.*cos(th).^2).*fTh;
            di = r./h.*cos(aol).*sin(inc).*fH_si;
            dO = r.*sin(aol)./(h).*fH_si;
            dw = -p./(h.*e).*cos(th).*fR +...
                r./(h.*e).*(2+e.*cos(th)).*sin(th).*fTh +...
                -r./h.*sin(aol).*cos(inc).*fH_si;
            dM = (p.*cos(th)-2*r.*e)./(n.*a.^2.*e).*fR +...
                -(p+r).*sin(th)./(n.*a.^2.*e).*fTh;
            
            dOe = reshape([da;de;di;dO;dw;dM + n],6*P.Con.nSats,1);
            % Equations of Motion
            dX = dOe;
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
        
        function dX = DynOeOsc3(P,t,X)  %#ok<INUSL>
            % rewrite with LPE
            OE  = reshape(X,6,P.Con.nSats);
            % Normalized Parameters
            mu = 1;
            Re = 1;
            % Element vectors (angles already in radians)
            a   = OE(1,:);
            e   = OE(2,:);
            inc = OE(3,:);
            w   = OE(5,:);
            Me  = OE(6,:);
            f  = pi/180*me2ta(Me*180/pi,e);
            % Secondary definitions
            p = a.*(1-e.^2);
            h = sqrt(mu*p);
            r = p./(1+e.*cos(f));
            k2 = -3*mu./r.^3.*P.Con.J2.*Re.^2;
            n = sqrt(mu./a.^3);
            eta = sqrt(1-e.^2);
            aol = f + w;
            % Test vs. GVE
%             MJ2 = -3*P.Con.mu./r.^4.*P.Con.J2.*P.Con.Re.^2;
%             % Forces
%             fR = MJ2/2.*(1-3*sin(inc).^2.*sin(aol).^2);
%             fTh = MJ2.*sin(inc).^2.*sin(aol).*cos(aol);
%             fH_si = MJ2.*cos(inc).*sin(aol); % sin(inc) taken out to avoid singularity
%             % Element Rates
%             da = 2*a.^2./h.*e.*sin(f).*fR +...
%                 2*a.^2./h.*(1 + e.*cos(f)).*fTh;
%             de = p./h.*sin(f).*fR +...
%                 r./h.*(e + 2*cos(f) + e.*cos(f).^2).*fTh;
%             di = r./h.*cos(aol).*sin(inc).*fH_si;
%             dO = r.*sin(aol)./(h).*fH_si;
%             dw = -p./(h.*e).*cos(f).*fR +...
%                 r./(h.*e).*(2+e.*cos(f)).*sin(f).*fTh +...
%                 -r./h.*sin(aol).*cos(inc).*fH_si;
%             dM = (p.*cos(f)-2*r.*e)./(n.*a.^2.*e).*fR +...
%                 -(p+r).*sin(f)./(n.*a.^2.*e).*fTh;

            dRde = ((e.^2.*cos(f)+cos(f)+2*e).*(3.*sin(inc).^2.*sin(aol).^2-1) ...
                - sin(f).*(2+e.*cos(f)).*(e.*sin(f).*(3*sin(inc).^2.*sin(aol).^2-1) ...
                - 2*sin(inc).^2.*sin(aol).*cos(aol).*(1+e.*cos(f))))./...
                (2*(1+e.*cos(f)).*(1-e.^2));
            % Element Rates
            da = k2.*h.*(e.*sin(f).*(1-3*sin(inc).^2.*sin(aol).^2) +...
                2*sin(inc).^2.*sin(aol).*cos(aol).*(1+e.*cos(f)))./(n.^2.*a.*p.*r);
            de = (1-e.^2).*da./(2.*a.*e) ...
                - k2.*eta.*sin(inc).^2.*sin(aol).*cos(aol)./(n.*a.^2.*e);
            di = k2.*sin(inc).*cos(inc).*cos(aol).*sin(aol)./(n.*a.^2.*eta);
            dO = k2.*cos(inc).*sin(aol).^2./(n.*a.^2.*eta);
            dw = k2.*eta.*dRde./(n.*a.^2.*e) - cos(inc).*dO;
            dM = k2.*(3.*sin(inc).^2.*sin(aol).^2-1)./(n.*a.^2) ...
                - k2.*eta.^2.*dRde./(n.*a.^2.*e);
            
            dOe = reshape([da;de;di;dO;dw;dM + n],6*P.Con.nSats,1);
            % Equations of Motion
            dX = dOe;
        end

        function dX = DynOePns(P,t,X)  %#ok<INUSL>
            % Dynamics for Nonsingular polar-nodal propagation
            mu = 1;
            J2 = P.Con.J2;
            Re = 1;
            PNS  = reshape(X,6,P.Con.nSats);
            % Element vectors (angles already in radians)
            r   = PNS(1,:); % radius
            th  = PNS(2,:); % mod. argument of latitude
            n   = PNS(3,:); % mod. right ascension
            R   = PNS(4,:); % Radial velocity
            TH  = PNS(5,:); % specific angular momentum
            N   = PNS(6,:); % mod. angular momentum in z
            
            % remove singularities n/N = nan iff N=0 i.e. orbit is equatorial
            nbyN = n./N;
            nbyN(isnan(nbyN)) = 0;
            k2 = J2*Re^2/2;
            
            % Keplerian Motion
            dKep = [R;
                    TH./r.^2;
                    zeros(1,P.Con.nSats);
                    TH.^2./r.^3-mu./r.^2;
                    zeros(1,P.Con.nSats);
                    zeros(1,P.Con.nSats)];
            den = (2*r.^3.*TH.^2);
            dJ2 = 3*mu*k2*[zeros(1,P.Con.nSats);
                 (N.^4-2*TH.*N.^2).*sin(th+nbyN).^2./(den.*TH);
                  -sin(th+nbyN).*((2*N.^3-4*TH.*N).*sin(th+nbyN) -cos(th+nbyN).*n.*(N.^2-4*TH))./(den);
                  -3*((N.^4-4*TH.*N.^2).*sin(th+nbyN).^2 +4*TH.^2/3)./(2*r.*den);
                  (N.^4-4*TH.*N.^2).*sin(th+nbyN).*cos(th+nbyN)./(den);
                  (N.^3-4*TH.*N).*sin(th+nbyN).*cos(th+nbyN)./(den)];
            dOe = reshape(dKep + dJ2,6*P.Con.nSats,1);
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
            
            J2 = P.Con.primary.J2;
            Re = P.Con.primary.Re;
            mu = P.Con.primary.mu;
            
            % handle elements vector
            a = X(1);
            e = X(2);
            i = X(3);
            aop = X(5);
            M = X(6);
            b = (1-sqrt(1-e^2))/e;
            
            
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
            
            
            C = [6*(3*sin(i)^2*cos(aop)^2-1)/(1-e^2)^2;
                (9*e^2*sin(i)^2*cos(aop)^2+3*sin(i)^2*(sin(aop)^2-cos(aop)^2)-3*e^2)/2/(1-e^2)^3.5;
                (3*e^3*sin(i)^2*cos(aop)^2+9*e*sin(i)^2*(sin(aop)^2-cos(aop)^2)-e^3)/4/(1-e^2)^4;
                (9*e^2*sin(i)^2*(sin(aop)^2-cos(aop)^2))/8/(1-e^2)^4.5;
                (3*e^3*sin(i)^2*(sin(aop)^2-cos(aop)^2))/16/(1-e^2)^5];
            
            dCdi_si = [36*cos(aop)^2*cos(i)/(1-e^2)^2; % all pre-divided by sin(i)
                3*cos(i)*(3*e^2*cos(aop)^2-2*cos(aop)^2+1)/(1-e^2)^3.5;
                3*e*cos(i)*(3+(e^2-6)*cos(aop)^2)/2/(1-e^2)^4;
                -9*e^2*cos(i)*(2*cos(aop)^2-1)/4/(1-e^2)^4.5;
                -3*e^3*cos(i)*(2*cos(aop)^2-1)/8/(1-e^2)^5];
            
            dCdo = [-36*sin(i)^2*sin(aop)*cos(aop)/(1-e^2)^2;
                3*(-3*e^2+2)*sin(i)^2*sin(aop)*cos(aop)/(1-e^2)^3.5;
                -3*e*(e^2-6)*sin(i)^2*sin(aop)*cos(aop)/2/(1-e^2)^4;
                9*e^2*sin(i)^2*sin(aop)*cos(aop)/2/(1-e^2)^4.5;
                3*e^3*sin(i)^2*sin(aop)*cos(aop)/4/(1-e^2)^5];
            
            dCdo_si = [-36*sin(i)*sin(aop)*cos(aop)/(1-e^2)^2;
                3*(-3*e^2+2)*sin(i)*sin(aop)*cos(aop)/(1-e^2)^3.5;
                -3*e*(e^2-6)*sin(i)*sin(aop)*cos(aop)/2/(1-e^2)^4;
                9*e^2*sin(i)*sin(aop)*cos(aop)/2/(1-e^2)^4.5;
                3*e^3*sin(i)*sin(aop)*cos(aop)/4/(1-e^2)^5];
            
            
            dCde = [24*e*(3*sin(i)^2*cos(aop)^2-1)/(1-e^2)^3;
                3*e*(15*e^2*sin(i)^2*cos(aop)^2 +7*sin(i)^2*sin(aop)^2 -sin(i)^2*cos(aop)^2 -5*e^2 -2)/2/(1-e^2)^4.5;
                ((15*e^4-54*e^2-9)*sin(i)^2*cos(aop)^2 + (63*e^2+9)*sin(i)^2*sin(aop)^2 -5*e^4 -3*e^2)/4/(1-e^2)^5;
                -9*e*sin(i)^2*(cos(aop)^2-sin(aop)^2)*(7*e^2+2)/8/(1-e^2)^5.5;
                -3*e^2*sin(i)^2*(cos(aop)^2-sin(aop)^2)*(7*e^2+3)/16/(1-e^2)^6];
            
            
            S = -[6*sin(i)^2*sin(aop)*cos(aop)/2/(1-e^2)^3;
                18*e*sin(i)^2*sin(aop)*cos(aop)/4/(1-e^2)^3.5;
                18*e^2*sin(i)^2*sin(aop)*cos(aop)/8/(1-e^2)^4;
                6*e^3*sin(i)^2*sin(aop)*cos(aop)/16/(1-e^2)^4.5];
            
            dSdi_si = -[6*cos(i)*sin(aop)*cos(aop)/(1-e^2)^3;
                9*e*cos(i)*sin(aop)*cos(aop)/(1-e^2)^3.5;
                9*e^2*cos(i)*sin(aop)*cos(aop)/2/(1-e^2)^4;
                3*e^3*cos(i)*sin(aop)*cos(aop)/4/(1-e^2)^4.5];
            
            dSdo = -[6*sin(i)^2*(2*cos(aop)^2-1)/2/(1-e^2)^3;
                18*e*sin(i)^2*(2*cos(aop)^2-1)/4/(1-e^2)^3.5;
                18*e^2*sin(i)^2*(2*cos(aop)^2-1)/8/(1-e^2)^4;
                6*e^3*sin(i)^2*(2*cos(aop)^2-1)/16/(1-e^2)^4.5];
            
            dSdo_si = -[6*sin(i)*(2*cos(aop)^2-1)/2/(1-e^2)^3;
                18*e*sin(i)*(2*cos(aop)^2-1)/4/(1-e^2)^3.5;
                18*e^2*sin(i)*(2*cos(aop)^2-1)/8/(1-e^2)^4;
                6*e^3*sin(i)*(2*cos(aop)^2-1)/16/(1-e^2)^4.5];
            
            dSde = [-18*e*sin(i)^2*sin(aop)*cos(aop)/(1-e^2)^4;
                -9*sin(i)^2*sin(aop)*cos(aop)*(6*e^2+1)/2/(1-e^2)^4.5;
                -9*e*sin(i)^2*sin(aop)*cos(aop)*(3*e^2+1)/2/(1-e^2)^5;
                -9*e^2*sin(i)^2*sin(aop)*cos(aop)*(2*e^2+1)/8/(1-e^2)^5.5];
            
            AkM = nan(5,kMax);
            Ak_eM = nan(5,kMax);
            Akde_eM = nan(5,kMax);
            
            BkM = nan(4,kMax);
            Bk_eM = nan(4,kMax);
            Bkde_eM = nan(4,kMax);
            
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
                
                dg2deXe = abs(m2+n+k-2).*b.^(abs(m2+n+k-2))/sqrt(1-e^2);
                dg3deXe = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^(abs(m3+n+k-3)) + ...
                    e/sqrt(1-e^2)*b.^(abs(m3+n+k-2)))/sqrt(1-e^2) + ...
                    e*b.^abs(m3+n+k-2)/(1-e^2)^(3/2);
                dg4deXe = abs(m4+n+k-2).*(3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2) +...
                    abs(m4+n+k-3).*(3*e/2/sqrt(1-e^2)*b.^abs(m4+n+k-3) +...
                    abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))/sqrt(1-e^2) +...
                    3/2*e*abs(m4+n+k-2)/(1-e^2)^(3/2).*b.^abs(m4+n+k-3) +...
                    3*e^2/(1-e^2)^2*b.^abs(m4+n+k-2);
                dg5deXe = abs(m5+n+k-2).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(m5+n+k-2) + ...
                    abs(m5+n+k-3).*(5*e^2/2/(1-e^2)*b.^abs(m5+n+k-3) + ...
                    abs(m5+n+k-4).*(e/sqrt(1-e^2)*b.^abs(m5+n+k-4) + ...
                    abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))/sqrt(1-e^2) + ...
                    abs(m5+n+k-3).*abs(m5+n+k-2)*e/(1-e^2)^(3/2).*b.^abs(m5+n+k-4) +...
                    5*abs(m5+n+k-2)*e^2/(1-e^2)^2.*b.^abs(m5+n+k-3) +...
                    15*e^3/2/(1-e^2)^(5/2)*b.^abs(m5+n+k-2);
                
                Jk = besselj(k,k*e);
                % Jk/e nonsingular
                Jk_e = 0.5*(besselj(k+1,k*e) + besselj(k-1,k*e));
                % dJkde/e nonsingular
                if k~=1
                    % use expression with no /e
                    dJkde_e = k^2/4/(k^2-1)*(2*Jk + (k+1)*besselj(k-2,k*e) - (k-1)*besselj(k+2,k*e));
                else
                    % elimination of /e not possible
                    dJkde_e = 0.5*k/e*(besselj(k-1,k*e) - besselj(k+1,k*e));
                end
                
                Jn = besselj(n,-k*e);
                % Jn/e nonsingular
                if n~=0
                    Jn_e = -k/2/n*(besselj(n+1,-k*e)+besselj(n-1,-k*e));
                else
                    Jn_e = Jn/e;
                end
                % dJnde/e nonsingular
                if abs(n)~=1
                    dJnde_e = k^2/4/(n^2-1)*(2*Jn + (n+1)*besselj(n-2,-k*e) - (n-1)*besselj(n+2,-k*e));
                else
                    dJnde_e = -0.5*k/e*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
                end
                % Jn/e^2 nonsingular
                if abs(n)>1
                    Jn_e2 = k^2/4/n/(n^2-1)*(2*n*Jn + (n+1)*besselj(n-2,-k*e) + (n-1)*besselj(n+2,-k*e));
                else
                    Jn_e2 = Jn/e^2;
                end
                
                Ak = [Jk; Jn*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                Ak_e = [Jk_e; Jn_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                Akde_e = [dJkde_e; dJnde_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] +...
                    Jn_e*[da2de.'*g2; da3de.'*g3; da4de.'*g4; da5de.'*g5] +...
                    Jn_e2*[a2.'*dg2deXe; a3.'*dg3deXe; a4.'*dg4deXe; a5.'*dg5deXe]];
                
                Bk = Jn*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                Bk_e = Jn_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                Bkde_e = dJnde_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] +...
                    Jn_e*[db1de.'*g2; db2de.'*g3; db3de.'*g4; db4de.'*g5] +...
                    Jn_e2*[b1.'*dg2deXe; b2.'*dg3deXe; b3.'*dg4deXe; b4.'*dg5deXe];
                
                n = 1;
                while n <= k + 5
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
                    
                    dg2deXe = abs(m2+n+k-2).*b.^(abs(m2+n+k-2))/sqrt(1-e^2);
                    dg3deXe = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^(abs(m3+n+k-3)) + ...
                        e/sqrt(1-e^2)*b.^(abs(m3+n+k-2)))/sqrt(1-e^2) + ...
                        e*b.^abs(m3+n+k-2)/(1-e^2)^(3/2);
                    dg4deXe = abs(m4+n+k-2).*(3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2) +...
                        abs(m4+n+k-3).*(3*e/2/sqrt(1-e^2)*b.^abs(m4+n+k-3) +...
                        abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))/sqrt(1-e^2) +...
                        3/2*e*abs(m4+n+k-2)/(1-e^2)^(3/2).*b.^abs(m4+n+k-3) +...
                        3*e^2/(1-e^2)^2*b.^abs(m4+n+k-2);
                    dg5deXe = abs(m5+n+k-2).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(m5+n+k-2) + ...
                        abs(m5+n+k-3).*(5*e^2/2/(1-e^2)*b.^abs(m5+n+k-3) + ...
                        abs(m5+n+k-4).*(e/sqrt(1-e^2)*b.^abs(m5+n+k-4) + ...
                        abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))/sqrt(1-e^2) + ...
                        abs(m5+n+k-3).*abs(m5+n+k-2)*e/(1-e^2)^(3/2).*b.^abs(m5+n+k-4) +...
                        5*abs(m5+n+k-2)*e^2/(1-e^2)^2.*b.^abs(m5+n+k-3) +...
                        15*e^3/2/(1-e^2)^(5/2)*b.^abs(m5+n+k-2);
                    
                    Jn = besselj(n,-k*e);
                    % Jn/e nonsingular
                    if n~=0
                        Jn_e = -k/2/n*(besselj(n+1,-k*e)+besselj(n-1,-k*e));
                    else
                        Jn_e = Jn/e;
                    end
                    % dJnde/e nonsingular
                    if abs(n)~=1
                        dJnde_e = k^2/4/(n^2-1)*(2*Jn + (n+1)*besselj(n-2,-k*e) - (n-1)*besselj(n+2,-k*e));
                    else
                        dJnde_e = -0.5*k/e*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
                    end
                    % Jn/e^2 nonsingular
                    if abs(n)>1
                        Jn_e2 = k^2/4/n/(n^2-1)*(2*n*Jn + (n+1)*besselj(n-2,-k*e) + (n-1)*besselj(n+2,-k*e));
                    else
                        Jn_e2 = Jn/e^2;
                    end
                    
                    dAk = [0; Jn*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                    dAk_e = [0; Jn_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                    dAkde_e = [0; dJnde_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] +...
                        Jn_e*[da2de.'*g2; da3de.'*g3; da4de.'*g4; da5de.'*g5] +...
                        Jn_e2*[a2.'*dg2deXe; a3.'*dg3deXe; a4.'*dg4deXe; a5.'*dg5deXe]];
                    
                    dBk = Jn*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBk_e = Jn_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBkde_e = dJnde_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] +...
                        Jn_e*[db1de.'*g2; db2de.'*g3; db3de.'*g4; db4de.'*g5] +...
                        Jn_e2*[b1.'*dg2deXe; b2.'*dg3deXe; b3.'*dg4deXe; b4.'*dg5deXe];
                    
                    
                    % negative n
                    n = -n;
                    g2 = b.^abs(m2+n+k-2);
                    g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
                    g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
                        3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
                        3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
                    g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
                        e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
                        5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
                        5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);
                    
                    dg2deXe = abs(m2+n+k-2).*b.^(abs(m2+n+k-2))/sqrt(1-e^2);
                    dg3deXe = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^(abs(m3+n+k-3)) + ...
                        e/sqrt(1-e^2)*b.^(abs(m3+n+k-2)))/sqrt(1-e^2) + ...
                        e*b.^abs(m3+n+k-2)/(1-e^2)^(3/2);
                    dg4deXe = abs(m4+n+k-2).*(3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2) +...
                        abs(m4+n+k-3).*(3*e/2/sqrt(1-e^2)*b.^abs(m4+n+k-3) +...
                        abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))/sqrt(1-e^2) +...
                        3/2*e*abs(m4+n+k-2)/(1-e^2)^(3/2).*b.^abs(m4+n+k-3) +...
                        3*e^2/(1-e^2)^2*b.^abs(m4+n+k-2);
                    dg5deXe = abs(m5+n+k-2).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(m5+n+k-2) + ...
                        abs(m5+n+k-3).*(5*e^2/2/(1-e^2)*b.^abs(m5+n+k-3) + ...
                        abs(m5+n+k-4).*(e/sqrt(1-e^2)*b.^abs(m5+n+k-4) + ...
                        abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))/sqrt(1-e^2) + ...
                        abs(m5+n+k-3).*abs(m5+n+k-2)*e/(1-e^2)^(3/2).*b.^abs(m5+n+k-4) +...
                        5*abs(m5+n+k-2)*e^2/(1-e^2)^2.*b.^abs(m5+n+k-3) +...
                        15*e^3/2/(1-e^2)^(5/2)*b.^abs(m5+n+k-2);
                    
                    Jn = besselj(n,-k*e);
                    % Jn/e nonsingular
                    if n~=0
                        Jn_e = -k/2/n*(besselj(n+1,-k*e)+besselj(n-1,-k*e));
                    else
                        Jn_e = Jn/e;
                    end
                    % dJnde/e nonsingular
                    if abs(n)~=1
                        dJnde_e = k^2/4/(n^2-1)*(2*Jn + (n+1)*besselj(n-2,-k*e) - (n-1)*besselj(n+2,-k*e));
                    else
                        dJnde_e = -0.5*k/e*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
                    end
                    % Jn/e^2 nonsingular
                    if abs(n)>1
                        Jn_e2 = k^2/4/n/(n^2-1)*(2*n*Jn + (n+1)*besselj(n-2,-k*e) + (n-1)*besselj(n+2,-k*e));
                    else
                        Jn_e2 = Jn/e^2;
                    end
                    
                    dAk = dAk + [0; Jn*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                    dAk_e = dAk_e + [0; Jn_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                    dAkde_e = dAkde_e + [0; dJnde_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] +...
                        Jn_e*[da2de.'*g2; da3de.'*g3; da4de.'*g4; da5de.'*g5] +...
                        Jn_e2*[a2.'*dg2deXe; a3.'*dg3deXe; a4.'*dg4deXe; a5.'*dg5deXe]];
                    
                    dBk = dBk + Jn*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBk_e = dBk_e + Jn_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBkde_e = dBkde_e + dJnde_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] +...
                        Jn_e*[db1de.'*g2; db2de.'*g3; db3de.'*g4; db4de.'*g5] +...
                        Jn_e2*[b1.'*dg2deXe; b2.'*dg3deXe; b3.'*dg4deXe; b4.'*dg5deXe];
                    
                    Ak = Ak + dAk;
                    Ak_e = Ak_e + dAk_e;
                    Akde_e = Akde_e + dAkde_e;
                    
                    Bk = Bk + dBk;
                    Bk_e = Bk_e + dBk_e;
                    Bkde_e = Bkde_e + dBkde_e;
                    
                    n = abs(n); % return n to positive value
                    n = n+1;
                end
                
                AkM(:,k) = Ak;
                Ak_eM(:,k) = Ak_e;
                Akde_eM(:,k) = Akde_e;
                
                BkM(:,k) = Bk;
                Bk_eM(:,k) = Bk_e;
                Bkde_eM(:,k) = Bkde_e;
                
                k = k+1;
            end
            
            
            n = sqrt(mu/a^3); % bad practice, should rename but whatever
            eta = sqrt(1-e^2);
            % constant potential values
            R = -n^2*Re^2/2; % common factor
            Rsma = -n*J2*Re^2/a;
            Recc = -n*eta*J2*Re^2/2/a^2;
            Rinc = -n*J2*Re^2*cos(i)/2/a^2/eta;
            Rran = -n*J2*Re^2/2/a^2/eta;
            Raop = -n*J2*Re^2/2/a^2;
            Rman = Raop;
            
            % Freq 0 elements without common factor
            ran0 = 3*cos(i)/eta^3;
            aop0 = -1.5*(5*cos(i)^2-1)/eta^4;
            man0 = -1.5*(3*cos(i)^2-1)/eta^3;
            
            % Calculate Spectrum of Elements
            k = 1:kMax;
            lpeSpec = nan(12,kMax+1);
            
            lpeSpec(1:2,:) = Rsma*[0, S.'*(BkM.*k);
                0, -C.'*(AkM.*k)];
            lpeSpec(3:4,:) = Recc*[0, eta*S.'*(Bk_eM.*k) - dCdo.'*Ak_eM;
                0, -eta*C.'*(Ak_eM.*k) - dSdo.'*Bk_eM];
            lpeSpec(5:6,:) = Rinc*[0, dCdo_si.'*AkM;
                0, dSdo_si.'*BkM];
            lpeSpec(7:8,:) = Rran*[ran0, dCdi_si.'*AkM;
                0, dSdi_si.'*BkM];
            lpeSpec(9:10,:) = Raop*[aop0, eta*(dCde.'*Ak_eM + C.'*Akde_eM) - cos(i)/eta*dCdi_si.'*AkM;
                0, eta*(dSde.'*Bk_eM + S.'*Bkde_eM) - cos(i)/eta*dSdi_si.'*BkM];
            lpeSpec(11:12,:) = Rman*[man0, -eta^2*(dCde.'*Ak_eM + C.'*Akde_eM) + 6*C.'*AkM;
                0, -eta^2*(dSde.'*Bk_eM + S.'*Bkde_eM) + 6*S.'*BkM];
            % ******* Output dX *******
            k = 0:kMax;
            Ck = cos(k.'*M);
            Sk = sin(k.'*M);
            
            dX = lpeSpec(1:2:11,:)*Ck + lpeSpec(2:2:12,:)*Sk;
            dX(6) = dX(6) + n;
            
        end
        
        function [freq0, lpeSpec] = DynOeFourier2Ord(P,t,icOsc,kMax) 
            %% Handle Input
            J2 = P.Con.primary.J2;
            Re = P.Con.primary.Re;
            mu = P.Con.primary.mu;
            
            nT = length(t);
            
            icM = osc2meSP(icOsc); % change to numerical mean
            icM(3:end) = icM(3:end)*pi/180;
            icOsc(3:end) = icOsc(3:end)*pi/180;
            % handle elements vector
            a = icM(1);
            e = icM(2);
            i = icM(3);
            
            nMo = sqrt(mu/a^3);
            p = a*(1-e^2);
            aop = icM(5) + (3/4*J2*(Re/p)^2*nMo*(5*cos(i)^2-1))*t;
            
            b = (1-sqrt(1-e^2))/e;
            
            %% constant vectors
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
            
            %% matrices by aop
            C = [6*(3*sin(i)^2*cos(aop).^2-1)/(1-e^2)^2;
                (9*e^2*sin(i)^2*cos(aop).^2+3*sin(i)^2*(sin(aop).^2-cos(aop).^2)-3*e^2)/2/(1-e^2)^3.5;
                (3*e^3*sin(i)^2*cos(aop).^2+9*e*sin(i)^2*(sin(aop).^2-cos(aop).^2)-e^3)/4/(1-e^2)^4;
                (9*e^2*sin(i)^2*(sin(aop).^2-cos(aop).^2))/8/(1-e^2)^4.5;
                (3*e^3*sin(i)^2*(sin(aop).^2-cos(aop).^2))/16/(1-e^2)^5];
            
            dCdi_si = [36*cos(aop).^2*cos(i)/(1-e^2)^2; % all pre-divided by sin(i)
                3*cos(i)*(3*e^2*cos(aop).^2-2*cos(aop).^2+1)/(1-e^2)^3.5;
                3*e*cos(i)*(3+(e^2-6)*cos(aop).^2)/2/(1-e^2)^4;
                -9*e^2*cos(i)*(2*cos(aop).^2-1)/4/(1-e^2)^4.5;
                -3*e^3*cos(i)*(2*cos(aop).^2-1)/8/(1-e^2)^5];
            
            dCdo = [-36*sin(i)^2*sin(aop).*cos(aop)/(1-e^2)^2;
                3*(-3*e^2+2)*sin(i)^2*sin(aop).*cos(aop)/(1-e^2)^3.5;
                -3*e*(e^2-6)*sin(i)^2*sin(aop).*cos(aop)/2/(1-e^2)^4;
                9*e^2*sin(i)^2*sin(aop).*cos(aop)/2/(1-e^2)^4.5;
                3*e^3*sin(i)^2*sin(aop).*cos(aop)/4/(1-e^2)^5];
            
            dCdo_si = [-36*sin(i)*sin(aop).*cos(aop)/(1-e^2)^2;
                3*(-3*e^2+2)*sin(i)*sin(aop).*cos(aop)/(1-e^2)^3.5;
                -3*e*(e^2-6)*sin(i)*sin(aop).*cos(aop)/2/(1-e^2)^4;
                9*e^2*sin(i)*sin(aop).*cos(aop)/2/(1-e^2)^4.5;
                3*e^3*sin(i)*sin(aop).*cos(aop)/4/(1-e^2)^5];
            
            dCde = [24*e*(3*sin(i)^2*cos(aop).^2-1)/(1-e^2)^3;
                3*e*(15*e^2*sin(i)^2*cos(aop).^2 +7*sin(i)^2*sin(aop).^2 -sin(i)^2*cos(aop).^2 -5*e^2 -2)/2/(1-e^2)^4.5;
                ((15*e^4-54*e^2-9)*sin(i)^2*cos(aop).^2 + (63*e^2+9)*sin(i)^2*sin(aop).^2 -5*e^4 -3*e^2)/4/(1-e^2)^5;
                -9*e*sin(i)^2*(cos(aop).^2-sin(aop).^2)*(7*e^2+2)/8/(1-e^2)^5.5;
                -3*e^2*sin(i)^2*(cos(aop).^2-sin(aop).^2)*(7*e^2+3)/16/(1-e^2)^6];
            
            S = -[6*sin(i)^2*sin(aop).*cos(aop)/2/(1-e^2)^3;
                18*e*sin(i)^2*sin(aop).*cos(aop)/4/(1-e^2)^3.5;
                18*e^2*sin(i)^2*sin(aop).*cos(aop)/8/(1-e^2)^4;
                6*e^3*sin(i)^2*sin(aop).*cos(aop)/16/(1-e^2)^4.5];
            
            dSdi_si = -[6*cos(i)*sin(aop).*cos(aop)/(1-e^2)^3;
                9*e*cos(i)*sin(aop).*cos(aop)/(1-e^2)^3.5;
                9*e^2*cos(i)*sin(aop).*cos(aop)/2/(1-e^2)^4;
                3*e^3*cos(i)*sin(aop).*cos(aop)/4/(1-e^2)^4.5];
            
            dSdo = -[6*sin(i)^2*(2*cos(aop).^2-1)/2/(1-e^2)^3;
                18*e*sin(i)^2*(2*cos(aop).^2-1)/4/(1-e^2)^3.5;
                18*e^2*sin(i)^2*(2*cos(aop).^2-1)/8/(1-e^2)^4;
                6*e^3*sin(i)^2*(2*cos(aop).^2-1)/16/(1-e^2)^4.5];
            
            dSdo_si = -[6*sin(i)*(2*cos(aop).^2-1)/2/(1-e^2)^3;
                18*e*sin(i)*(2*cos(aop).^2-1)/4/(1-e^2)^3.5;
                18*e^2*sin(i)*(2*cos(aop).^2-1)/8/(1-e^2)^4;
                6*e^3*sin(i)*(2*cos(aop).^2-1)/16/(1-e^2)^4.5];
            
            dSde = [-18*e*sin(i)^2*sin(aop).*cos(aop)/(1-e^2)^4;
                -9*sin(i)^2*sin(aop).*cos(aop)*(6*e^2+1)/2/(1-e^2)^4.5;
                -9*e*sin(i)^2*sin(aop).*cos(aop)*(3*e^2+1)/2/(1-e^2)^5;
                -9*e^2*sin(i)^2*sin(aop).*cos(aop)*(2*e^2+1)/8/(1-e^2)^5.5];
            
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
            AkM = nan(5,kMax);
            Ak_eM = nan(5,kMax);
            Akde_eM = nan(5,kMax);
            
            BkM = nan(4,kMax);
            Bk_eM = nan(4,kMax);
            Bkde_eM = nan(4,kMax);
            
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
                
                dg2deXe = abs(m2+n+k-2).*b.^(abs(m2+n+k-2))/sqrt(1-e^2);
                dg3deXe = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^(abs(m3+n+k-3)) + ...
                    e/sqrt(1-e^2)*b.^(abs(m3+n+k-2)))/sqrt(1-e^2) + ...
                    e*b.^abs(m3+n+k-2)/(1-e^2)^(3/2);
                dg4deXe = abs(m4+n+k-2).*(3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2) +...
                    abs(m4+n+k-3).*(3*e/2/sqrt(1-e^2)*b.^abs(m4+n+k-3) +...
                    abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))/sqrt(1-e^2) +...
                    3/2*e*abs(m4+n+k-2)/(1-e^2)^(3/2).*b.^abs(m4+n+k-3) +...
                    3*e^2/(1-e^2)^2*b.^abs(m4+n+k-2);
                dg5deXe = abs(m5+n+k-2).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(m5+n+k-2) + ...
                    abs(m5+n+k-3).*(5*e^2/2/(1-e^2)*b.^abs(m5+n+k-3) + ...
                    abs(m5+n+k-4).*(e/sqrt(1-e^2)*b.^abs(m5+n+k-4) + ...
                    abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))/sqrt(1-e^2) + ...
                    abs(m5+n+k-3).*abs(m5+n+k-2)*e/(1-e^2)^(3/2).*b.^abs(m5+n+k-4) +...
                    5*abs(m5+n+k-2)*e^2/(1-e^2)^2.*b.^abs(m5+n+k-3) +...
                    15*e^3/2/(1-e^2)^(5/2)*b.^abs(m5+n+k-2);
                
                Jk = besselj(k,k*e);
                % Jk/e nonsingular
                Jk_e = 0.5*(besselj(k+1,k*e) + besselj(k-1,k*e));
                % dJkde/e nonsingular
                if k~=1
                    % use expression with no /e
                    dJkde_e = k^2/4/(k^2-1)*(2*Jk + (k+1)*besselj(k-2,k*e) - (k-1)*besselj(k+2,k*e));
                else
                    % elimination of /e not possible
                    dJkde_e = 0.5*k/e*(besselj(k-1,k*e) - besselj(k+1,k*e));
                end
                
                Jn = besselj(n,-k*e);
                % Jn/e nonsingular
                if n~=0
                    Jn_e = -k/2/n*(besselj(n+1,-k*e)+besselj(n-1,-k*e));
                else
                    Jn_e = Jn/e;
                end
                % dJnde/e nonsingular
                if abs(n)~=1
                    dJnde_e = k^2/4/(n^2-1)*(2*Jn + (n+1)*besselj(n-2,-k*e) - (n-1)*besselj(n+2,-k*e));
                else
                    dJnde_e = -0.5*k/e*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
                end
                % Jn/e^2 nonsingular
                if abs(n)>1
                    Jn_e2 = k^2/4/n/(n^2-1)*(2*n*Jn + (n+1)*besselj(n-2,-k*e) + (n-1)*besselj(n+2,-k*e));
                else
                    Jn_e2 = Jn/e^2;
                end
                
                Ak = [Jk; Jn*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                Ak_e = [Jk_e; Jn_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                Akde_e = [dJkde_e; dJnde_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] +...
                    Jn_e*[da2de.'*g2; da3de.'*g3; da4de.'*g4; da5de.'*g5] +...
                    Jn_e2*[a2.'*dg2deXe; a3.'*dg3deXe; a4.'*dg4deXe; a5.'*dg5deXe]];
                
                Bk = Jn*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                Bk_e = Jn_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                Bkde_e = dJnde_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] +...
                    Jn_e*[db1de.'*g2; db2de.'*g3; db3de.'*g4; db4de.'*g5] +...
                    Jn_e2*[b1.'*dg2deXe; b2.'*dg3deXe; b3.'*dg4deXe; b4.'*dg5deXe];
                
                n = 1;
                while n <= k + 5
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
                    
                    dg2deXe = abs(m2+n+k-2).*b.^(abs(m2+n+k-2))/sqrt(1-e^2);
                    dg3deXe = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^(abs(m3+n+k-3)) + ...
                        e/sqrt(1-e^2)*b.^(abs(m3+n+k-2)))/sqrt(1-e^2) + ...
                        e*b.^abs(m3+n+k-2)/(1-e^2)^(3/2);
                    dg4deXe = abs(m4+n+k-2).*(3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2) +...
                        abs(m4+n+k-3).*(3*e/2/sqrt(1-e^2)*b.^abs(m4+n+k-3) +...
                        abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))/sqrt(1-e^2) +...
                        3/2*e*abs(m4+n+k-2)/(1-e^2)^(3/2).*b.^abs(m4+n+k-3) +...
                        3*e^2/(1-e^2)^2*b.^abs(m4+n+k-2);
                    dg5deXe = abs(m5+n+k-2).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(m5+n+k-2) + ...
                        abs(m5+n+k-3).*(5*e^2/2/(1-e^2)*b.^abs(m5+n+k-3) + ...
                        abs(m5+n+k-4).*(e/sqrt(1-e^2)*b.^abs(m5+n+k-4) + ...
                        abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))/sqrt(1-e^2) + ...
                        abs(m5+n+k-3).*abs(m5+n+k-2)*e/(1-e^2)^(3/2).*b.^abs(m5+n+k-4) +...
                        5*abs(m5+n+k-2)*e^2/(1-e^2)^2.*b.^abs(m5+n+k-3) +...
                        15*e^3/2/(1-e^2)^(5/2)*b.^abs(m5+n+k-2);
                    
                    Jn = besselj(n,-k*e);
                    % Jn/e nonsingular
                    if n~=0
                        Jn_e = -k/2/n*(besselj(n+1,-k*e)+besselj(n-1,-k*e));
                    else
                        Jn_e = Jn/e;
                    end
                    % dJnde/e nonsingular
                    if abs(n)~=1
                        dJnde_e = k^2/4/(n^2-1)*(2*Jn + (n+1)*besselj(n-2,-k*e) - (n-1)*besselj(n+2,-k*e));
                    else
                        dJnde_e = -0.5*k/e*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
                    end
                    % Jn/e^2 nonsingular
                    if abs(n)>1
                        Jn_e2 = k^2/4/n/(n^2-1)*(2*n*Jn + (n+1)*besselj(n-2,-k*e) + (n-1)*besselj(n+2,-k*e));
                    else
                        Jn_e2 = Jn/e^2;
                    end
                    
                    dAk = [0; Jn*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                    dAk_e = [0; Jn_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                    dAkde_e = [0; dJnde_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] +...
                        Jn_e*[da2de.'*g2; da3de.'*g3; da4de.'*g4; da5de.'*g5] +...
                        Jn_e2*[a2.'*dg2deXe; a3.'*dg3deXe; a4.'*dg4deXe; a5.'*dg5deXe]];
                    
                    dBk = Jn*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBk_e = Jn_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBkde_e = dJnde_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] +...
                        Jn_e*[db1de.'*g2; db2de.'*g3; db3de.'*g4; db4de.'*g5] +...
                        Jn_e2*[b1.'*dg2deXe; b2.'*dg3deXe; b3.'*dg4deXe; b4.'*dg5deXe];
                    
                    
                    % negative n
                    n = -n;
                    g2 = b.^abs(m2+n+k-2);
                    g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
                    g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
                        3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
                        3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
                    g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
                        e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
                        5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
                        5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);
                    
                    dg2deXe = abs(m2+n+k-2).*b.^(abs(m2+n+k-2))/sqrt(1-e^2);
                    dg3deXe = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^(abs(m3+n+k-3)) + ...
                        e/sqrt(1-e^2)*b.^(abs(m3+n+k-2)))/sqrt(1-e^2) + ...
                        e*b.^abs(m3+n+k-2)/(1-e^2)^(3/2);
                    dg4deXe = abs(m4+n+k-2).*(3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2) +...
                        abs(m4+n+k-3).*(3*e/2/sqrt(1-e^2)*b.^abs(m4+n+k-3) +...
                        abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))/sqrt(1-e^2) +...
                        3/2*e*abs(m4+n+k-2)/(1-e^2)^(3/2).*b.^abs(m4+n+k-3) +...
                        3*e^2/(1-e^2)^2*b.^abs(m4+n+k-2);
                    dg5deXe = abs(m5+n+k-2).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(m5+n+k-2) + ...
                        abs(m5+n+k-3).*(5*e^2/2/(1-e^2)*b.^abs(m5+n+k-3) + ...
                        abs(m5+n+k-4).*(e/sqrt(1-e^2)*b.^abs(m5+n+k-4) + ...
                        abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))/sqrt(1-e^2) + ...
                        abs(m5+n+k-3).*abs(m5+n+k-2)*e/(1-e^2)^(3/2).*b.^abs(m5+n+k-4) +...
                        5*abs(m5+n+k-2)*e^2/(1-e^2)^2.*b.^abs(m5+n+k-3) +...
                        15*e^3/2/(1-e^2)^(5/2)*b.^abs(m5+n+k-2);
                    
                    Jn = besselj(n,-k*e);
                    % Jn/e nonsingular
                    if n~=0
                        Jn_e = -k/2/n*(besselj(n+1,-k*e)+besselj(n-1,-k*e));
                    else
                        Jn_e = Jn/e;
                    end
                    % dJnde/e nonsingular
                    if abs(n)~=1
                        dJnde_e = k^2/4/(n^2-1)*(2*Jn + (n+1)*besselj(n-2,-k*e) - (n-1)*besselj(n+2,-k*e));
                    else
                        dJnde_e = -0.5*k/e*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
                    end
                    % Jn/e^2 nonsingular
                    if abs(n)>1
                        Jn_e2 = k^2/4/n/(n^2-1)*(2*n*Jn + (n+1)*besselj(n-2,-k*e) + (n-1)*besselj(n+2,-k*e));
                    else
                        Jn_e2 = Jn/e^2;
                    end
                    
                    dAk = dAk + [0; Jn*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                    dAk_e = dAk_e + [0; Jn_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                    dAkde_e = dAkde_e + [0; dJnde_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] +...
                        Jn_e*[da2de.'*g2; da3de.'*g3; da4de.'*g4; da5de.'*g5] +...
                        Jn_e2*[a2.'*dg2deXe; a3.'*dg3deXe; a4.'*dg4deXe; a5.'*dg5deXe]];
                    
                    dBk = dBk + Jn*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBk_e = dBk_e + Jn_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBkde_e = dBkde_e + dJnde_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] +...
                        Jn_e*[db1de.'*g2; db2de.'*g3; db3de.'*g4; db4de.'*g5] +...
                        Jn_e2*[b1.'*dg2deXe; b2.'*dg3deXe; b3.'*dg4deXe; b4.'*dg5deXe];
                    
                    Ak = Ak + dAk;
                    Ak_e = Ak_e + dAk_e;
                    Akde_e = Akde_e + dAkde_e;
                    
                    Bk = Bk + dBk;
                    Bk_e = Bk_e + dBk_e;
                    Bkde_e = Bkde_e + dBkde_e;
                    
                    n = abs(n); % return n to positive value
                    n = n+1;
                end
                
                AkM(:,k) = Ak;
                Ak_eM(:,k) = Ak_e;
                Akde_eM(:,k) = Akde_e;
                
                BkM(:,k) = Bk;
                Bk_eM(:,k) = Bk_e;
                Bkde_eM(:,k) = Bkde_e;
                
                k = k+1;
            end
            %% Define Final constants
            eta = sqrt(1-e^2);
            % constant potential values
            R = -nMo^2*Re^2/2; % common factor
            Rsma = -nMo*J2*Re^2/a;
            Recc = -nMo*eta*J2*Re^2/2/a^2;
            Rinc = -nMo*J2*Re^2*cos(i)/2/a^2/eta;
            Rran = -nMo*J2*Re^2/2/a^2/eta;
            Raop = -nMo*J2*Re^2/2/a^2;
            Rman = Raop;
            
            % Freq 0 elements without common factor
            ran0 = 3*cos(i)/eta^3;
            aop0 = -1.5*(5*cos(i)^2-1)/eta^4;
            man0 = -1.5*(3*cos(i)^2-1)/eta^3;
            
            % Calculate Spectrum of Elements
            k = 1:kMax;
            lpeSpec = nan(6*2*kMax,nT);
            
            % first order stuff
            
%             lpeSpec(1:2,:) = Rsma*[0, S.'*(BkM.*k);
%                 0, -C.'*(AkM.*k)];
%             lpeSpec(3:4,:) = Recc*[0, eta*S.'*(Bk_eM.*k) - dCdo.'*Ak_eM;
%                 0, -eta*C.'*(Ak_eM.*k) - dSdo.'*Bk_eM];
%             lpeSpec(5:6,:) = Rinc*[0, dCdo_si.'*AkM;
%                 0, dSdo_si.'*BkM];
%             lpeSpec(7:8,:) = Rran*[ran0, dCdi_si.'*AkM;
%                 0, dSdi_si.'*BkM];
%             lpeSpec(9:10,:) = Raop*[aop0, eta*(dCde.'*Ak_eM + C.'*Akde_eM) - cos(i)/eta*dCdi_si.'*AkM;
%                 0, eta*(dSde.'*Bk_eM + S.'*Bkde_eM) - cos(i)/eta*dSdi_si.'*BkM];
%             lpeSpec(11:12,:) = Rman*[man0, -eta^2*(dCde.'*Ak_eM + C.'*Akde_eM) + 6*C.'*AkM;
%                 0, -eta^2*(dSde.'*Bk_eM + S.'*Bkde_eM) + 6*S.'*BkM];
%             
            %% Integrate 2nd Order components
            % dAop = 0*0.75*J2*Re^2/p^2*(5*cos(i)^2-1); %  O(J2^2)
            
            dSmaSpec = Rsma*[S.'*(BkM.*k),...
                -C.'*(AkM.*k)].';
            
            dEccSpec = Recc*[eta*S.'*(Bk_eM.*k) - dCdo.'*Ak_eM,...
                -eta*C.'*(Ak_eM.*k) - dSdo.'*Bk_eM].';
            
            dIncSpec = Rinc*[dCdo_si.'*AkM,...
                dSdo_si.'*BkM].';
            
            dRanSpec = Rran*[dCdi_si.'*AkM,...
                dSdi_si.'*BkM].';
            
            
            dAopSpec = Raop*[eta*(dCde.'*Ak_eM + C.'*Akde_eM) - cos(i)/eta*dCdi_si.'*AkM,...
                eta*(dSde.'*Bk_eM + S.'*Bkde_eM) - cos(i)/eta*dSdi_si.'*BkM].';
            
            dManSpec = Rman*[-eta^2*(dCde.'*Ak_eM + C.'*Akde_eM) + 6*C.'*AkM,...
                -eta^2*(dSde.'*Bk_eM + S.'*Bkde_eM) + 6*S.'*BkM].';
            
            lpeSpec = [dSmaSpec;
                dEccSpec;
                dIncSpec;
                dRanSpec;
                dAopSpec;
                dManSpec];
            
            %% Zero Frequency Components
            freq0 = [0; 0; 0; Rran*ran0; Raop*aop0; Rman*man0];
            
            
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