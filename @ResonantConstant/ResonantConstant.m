classdef ResonantConstant < OrbitControl
    % Constant magnitude resonant controller for a, e, i
    % Based on work with Bharath
    properties
        acc % thrust acceleration magnitude
        targetOe % Target orbital element vector - [a;e;i]
        tol % Tolerance to shut off controller
    end

    methods
        function RCC = ResonantConstant(TargetOe,ThrustMagnitude,TargetTolerance,primary)
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
            RCC = RCC@OrbitControl(primary);
            
            %%%% Post Initialization %%%%
            % property assignment
            RCC.acc = ThrustMagnitude;
            if size(TargetOe,1) ~= 3
                TargetOe = TargetOe.';
            end
            RCC.targetOe = TargetOe;
            RCC.tol = TargetTolerance;
        end

        function fRsw = ControlRSW(RCC,t,oeOsc)
            nT = size(oeOsc,2);
            % output control acceleration in RSW frame
            % Normalized Units
            % oeOsc - 6 x nT, normalized units
            % Scaling
            mu = 1;
            Re = 1;
            dScale = RCC.primary.Re;
            tScale = sqrt(RCC.primary.Re^3/RCC.primary.mu);
            fScale = dScale/tScale^2;
            % Handle OE - controller uses mean?
            oeOsc(3:end,:) = oeOsc(3:end,:)*180/pi;
            oeM = osc2me(oeOsc,RCC.primary.J2,Re);
            sma = oeM(1,:);
            ecc = oeM(2,:);
            inc = oeM(3,:)*pi/180;
            aop = oeM(5,:)*pi/180;
            tan = me2ta(oeM(6,:),oeM(2,:))*pi/180;
            

            dOe = [sma;ecc;inc] - RCC.targetOe./[dScale,1,180/pi].'; % mean element error 3 x nT
            indTol = (vecnorm(dOe,2,1)/norm(RCC.targetOe./[dScale,1,180/pi].')<RCC.tol);
            if nT == 1 && indTol
                fRsw = zeros(3,1);
            else
            % Controller from Bharath
            if(any(indTol))
            fRsw(:,indTol) = zeros(3,1);
            end
            fRsw = (RCC.acc/fScale)*fRsw./vecnorm(fRsw,2,1);
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
