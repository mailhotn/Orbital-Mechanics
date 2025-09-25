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

        function [Time, X] = PropOeMeanFast(P,T)
            % Refactored - legacy
            [Time, X] = PropOeMean(P,T);
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
            % Use OeOsc3 probably
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
            % Numerically propagate LPE for oscullating elements
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
            % Very stupid - only for demonstrating that Fourier series is
            % correct.
            % Currently only works for one satellite
            opts = odeset('reltol',P.relTol,'abstol',P.absTol);
            IC = reshape(P.Con.InitialOeOsc,[6*P.Con.nSats,1]);
            IC(3:end) = IC(3:end)*pi/180;
            [Time, X] = ode45(@(T,X) P.DynOeFourier(T,X,kMax),T,IC,opts);
            X(:,3:end) = X(:,3:end)*180/pi;
            %             X = X;
            %             X(3:end,:) = wrapTo360(X(3:end,:));
        end

        function [Time, X] = PropOeFourierNoMFix(P,T,kMax)
            % Retest if sequential M is good. It's less rigorous?
            % Propagate for time T using Fourier LPE
            % Assume constant coefficients

            icVec = reshape(P.Con.InitialOeOsc,[6*P.Con.nSats,1]); % ic of all sats
            X = nan(6*P.Con.nSats,length(T));
            % for satellites - inefficient, could maybe get spectrum of
            % multiple satellites, not sure if worth the effort though
            for iSat = 1:P.Con.nSats
                icOsc = icVec((1:6)+6*(iSat-1));
                [freq0,lpeSpec] = P.DynOeFourier([],icOsc,kMax);
                icM = osc2meNum(icOsc);
                icM(3:end) = icM(3:end)*pi/180;
                icOsc(3:end) = icOsc(3:end)*pi/180;

                smaM = icM(1);

                nM = sqrt(P.Con.primary.mu/smaM^3);
                M = nM*T+icOsc(6);
                k = 1:kMax;

                Xi = nan(6,length(T));

                % Calculate all elements
                Sk = sin(k.'*M)./k.'/nM;
                Ck = -cos(k.'*M)./k.'/nM;
                Ak = lpeSpec(1:2:11,:);
                Bk = lpeSpec(2:2:12,:);
                fourIntSolAll = Ak*Sk + Bk*Ck;
                % Assign values
                Xi = icOsc + freq0*T + fourIntSolAll - fourIntSolAll(:,1);
                Xi(6,:) = Xi(6,:) + M - icOsc(6);


                Xi(3:5,:) = wrapTo360(Xi(3:5,:)*180/pi);
                Xi(6,:) = Xi(6,:)*180/pi;
                % finished Xi, assign then iterate
                X((1:6)+6*(iSat-1),:) = Xi;
            end
            X = X.';
            Time = T;
        end

        function [Time, X] = PropOeFourier(P,T,kMax)
            % Propagate for time T using Fourier LPE
            % Assume constant coefficients

            %% Handle Initial conditions
            icVec = reshape(P.Con.InitialOeOsc,[6*P.Con.nSats,1]); % ic of all sats
            X = nan(6*P.Con.nSats,length(T));
            % for satellites - inefficient, could maybe get spectrum of
            % multiple satellites, not sure if worth the effort though
            for iSat = 1:P.Con.nSats
                icOsc = icVec((1:6)+6*(iSat-1));
                [freq0,lpeSpec] = P.DynOeFourier([],icOsc,kMax);
                icM = osc2meNum(icOsc);
                icM(3:end) = icM(3:end)*pi/180;
                icOsc(3:end) = icOsc(3:end)*pi/180;

                smaM = icM(1);

                nM = sqrt(P.Con.primary.mu/smaM^3);
                %             n = sqrt(P.Con.primary.mu/a^3)*(1+3*g2*(1-1.5*sin(i)^2)/eta^3); % kozai Fix
                %                         n = n + lpeSpec(11,1); % <-------------------  Work on this
                M = nM*T+icOsc(6);
                k = 1:kMax;

                Xi = nan(6,length(T));

                % Fix M - add first Fourier variations
                Sk = sin(k.'*M)./k.'/nM;
                Ck = -cos(k.'*M)./k.'/nM;
                AkM = lpeSpec(11,:);
                BkM = lpeSpec(12,:);
                fourIntSolM = AkM*Sk + BkM*Ck;
                M2 = freq0(6)*T + fourIntSolM - fourIntSolM(1) + M;

                % Calculate all other elements
                Sk = sin(k.'*M2)./k.'/nM;
                Ck = -cos(k.'*M2)./k.'/nM;
                Ak = lpeSpec(1:2:11,:);
                Bk = lpeSpec(2:2:12,:);
                fourIntSolAll = Ak*Sk + Bk*Ck;

                Xi = icOsc + freq0*T + fourIntSolAll - fourIntSolAll(:,1);
                Xi(6,:) = Xi(6,:) + M - icOsc(6);


                Xi(3:5,:) = wrapTo360(Xi(3:5,:)*180/pi);
                Xi(6,:) = Xi(6,:)*180/pi;
                % finished Xi, assign then iterate
                X((1:6)+6*(iSat-1),:) = Xi;
            end
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
                icM = osc2meNum(icOsc); % changed to numerical mean
                icOsc(3:end) = icOsc(3:end)*pi/180;
                icM(3:end) = icM(3:end)*pi/180;

                smaM = icM(1);

                % Continue with the rest
                nM = sqrt(P.Con.primary.mu/smaM^3);

                M = nM*T+icOsc(6);
                k = 1:kMax;
                Xi = nan(6,length(T)); % Xi is result for current iteration

                % % Initial M - Possibly redundant
                % trigMat = [sin(k.'*M(1))./k.'/nM;-cos(k.'*M(1))./k.'/nM]; % 2kx1 - Sk, Ck for t=0
                % InitM = sum(lpeSpec((10*kMax+1):end,1).*trigMat);

                % Fix M - from linear M to variations based off Fourier
                Sk = sin(k.'*M)./k.'/nM; % kxnT
                Ck = -cos(k.'*M)./k.'/nM; % kxnT
                AkM = lpeSpec((10*kMax+1):11*kMax,:);
                BkM = lpeSpec((11*kMax+1):12*kMax,:);
                fourIntSolM = sum(AkM.*Sk + BkM.*Ck);
                M2 = freq0(6)*T + fourIntSolM - fourIntSolM(1) + M;

                % % Initial time - With corrected M
                % trigMat = repmat([sin(k.'*M2(1))./k.'/nM;-cos(k.'*M2(1))./k.'/nM],6,1);
                % trigsum1 = lpeSpec(:,1).*trigMat;
                % % CHANGE TO MATRIX
                % % MULTIPLICATION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                % InitVal = [sum(trigsum1(1:2*kMax)); sum(trigsum1((2*kMax+1):4*kMax));...
                %     sum(trigsum1((4*kMax+1):6*kMax)); sum(trigsum1((6*kMax+1):8*kMax));...
                %     sum(trigsum1((8*kMax+1):10*kMax)); sum(trigsum1((10*kMax+1):12*kMax))];


                % Calculate all elements
                Sk = repmat(sin(k.'*M2)./k.'/nM,6,1);
                Ck = repmat(-cos(k.'*M2)./k.'/nM,6,1);
                iAk = (1:(6*kMax)) + reshape(repmat((0:5)*kMax,kMax,1),1,kMax*6); % Indices for Ak in Spectrum
                iBk = (1:(6*kMax)) + reshape(repmat((1:6)*kMax,kMax,1),1,kMax*6);
                Ak = lpeSpec(iAk,:);
                Bk = lpeSpec(iBk,:);
                sumMat = blkdiag(ones(1,kMax),ones(1,kMax),ones(1,kMax),ones(1,kMax),ones(1,kMax),ones(1,kMax));
                fourIntSolAll = sumMat*(Ak.*Sk + Bk.*Ck); % Solution to integral of LPE in Fourier Form

                Xi = icOsc + freq0*T + fourIntSolAll - fourIntSolAll(:,1);
                %             X(6,:) = M2;
                Xi(6,:) = Xi(6,:) + M - icOsc(6); % subtracted M(0)
                % why M instead of M2? because M2 would be adding variations twice


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
            amoP = sqrt((1-ecc^2)*sma*mu);            % Theta - angular momentum
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

        dX = DynEciJ2(P,t,X)

        dX = DynOeOsc(P,t,X) 

        dX = DynOeOsc2(P,t,X) 

        dX = DynOeOsc3(P,t,X)

        dX = DynOePns(P,t,X) 

        dX = DynEciJ3(P,t,X)
        
        [freq0,lpeSpec] = DynOeFourier(P,t,icOsc,kMax)

        [freq0, lpeSpec] = DynOeFourier2Ord(P,t,icOsc,kMax)
    
    end

    methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(P)  %#ok<MANU>
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Propagation Definitions';
            propgroups(1).PropertyList = {'relTol','absTol','Con'};
        end
    end
end