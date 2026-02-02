classdef ZhangInverse < OrbitControl
    % Constant magnitude inverse dynamics controller for a, e, i
    % H. Zhang and P. Gurfil, “Nanosatellite Cluster Keeping Under Thrust
    % Uncertainties," Journal of Guidance, Control, and Dynamics, Sep. 2014, doi: 10.2514/1.G000554.
    % Can compare to Constant mag resonant
    properties
        acc % thrust acceleration magnitude
        targetOe % Target orbital element vector - [a;e;i]
        tol % Tolerance to shut off controller
    end

    methods
        function ZI = ZhangInverse(TargetOe,ThrustMagnitude,TargetTolerance,primary)
            %%%% Pre Initialization %%%%
            switch nargin
                case 0 % default tangential
                    error('Need Target')
                case 1 % default thrust mag
                    ThrustMagnitude = 0.00001; % 0.001 m/s^2
                    TargetTolerance = 1e-5;
                    primary = Earth;
                case 2 % arbitrary thrust mag
                    TargetTolerance = 1e-5;
                    primary = Earth;
                case 3
                    primary = Earth;    
                case 4

                otherwise
                    error('Wrong number of input arguments')
            end
            %%%% Object Initialization %%%%
            % Call superclass constructor before accessing object
            ZI = ZI@OrbitControl(primary);
            
            %%%% Post Initialization %%%%
            % property assignment
            ZI.acc = ThrustMagnitude;
            if size(TargetOe,1) ~= 3
                TargetOe = TargetOe.';
            end
            ZI.targetOe = TargetOe;
            ZI.tol = TargetTolerance;
        end

        function fRsw = ControlRSW(ZI,t,oeOsc)
            nT = size(oeOsc,2);
            % output control acceleration in RSW frame
            % Normalized Units
            % oeOsc - 6 x nT, normalized units
            % Scaling
            mu = 1;
            Re = 1;
            dScale = ZI.primary.Re;
            tScale = sqrt(ZI.primary.Re^3/ZI.primary.mu);
            fScale = dScale/tScale^2;
            % Handle OE - controller uses mean
            oeOsc(3:end,:) = oeOsc(3:end,:)*180/pi;
            oeM = osc2me(oeOsc,ZI.primary.J2,Re);
            sma = oeM(1,:);
            ecc = oeM(2,:);
            inc = oeM(3,:)*pi/180;
            aop = oeM(5,:)*pi/180;
            tan = me2ta(oeM(6,:),oeM(2,:))*pi/180;
            

            dOe = [sma;ecc;inc] - ZI.targetOe./[dScale,1,180/pi].'; % mean element error 3 x nT
            indTol = (vecnorm(dOe,2,1)/norm(ZI.targetOe./[dScale,1,180/pi].')<ZI.tol);
            if nT == 1 && indTol
                fRsw = zeros(3,1);
            else
            dOe = reshape(dOe,3*nT,1); % 3*nT x 1 for multiplication
            % Get controller matrix
            eta = sqrt(1-ecc.^2);
            B11 = 1/sqrt(mu)*2*sma.^(3/2).*ecc.*sin(tan)./eta;
            B12 = 1/sqrt(mu)*2*sma.^(3/2).*(ecc.*cos(tan)+1)./eta;
            % B13 = zeros(1,nT);
            B21 = 1/sqrt(mu)*sqrt(sma).*eta.*sin(tan);
            B22 = 1/sqrt(mu)*sqrt(sma).*eta.*(cos(tan).*(ecc.*cos(tan)+2)+ecc)./(ecc.*cos(tan)+1);
            % B23 = zeros(1,nT);
            % B31 = zeros(1,nT);
            % B32 = zeros(1,nT);
            B33 = 1/sqrt(mu)*sqrt(sma).*eta.*cos(tan+aop)./(ecc.*cos(tan)+1);
            
            D1 = diag(reshape([B11;B22;B33],1,nT*3));
            D2 = reshape([B12;zeros(2,nT)],1,nT*3);
            D2(end) = [];
            D2 = diag(D2,1);
            D3 = reshape([B21;zeros(2,nT)],1,nT*3);
            D3(end) = [];
            D3 = diag(D3,-1);
            B = D1+D2+D3;

            % Calculate control acceleration
            fRsw = -B\dOe; % 3*nT x 1
            fRsw = reshape(fRsw,3,nT);
            if(any(indTol))
            fRsw(:,indTol) = zeros(3,1);
            end
            fRsw = (ZI.acc/fScale)*fRsw./vecnorm(fRsw,2,1);
            end
            
        end
    end
    



    methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(Sat) %#ok
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Thrust Magnitude';
            propgroups(1).PropertyList = {'acc'};
            propgroups(2).Title = 'Thrust Direction';
            propgroups(2).PropertyList = {'dir'};

        end
    end
end
