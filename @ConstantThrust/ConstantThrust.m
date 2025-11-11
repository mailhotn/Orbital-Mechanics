classdef ConstantThrust < OrbitControl
    % Constant thrust in rsw frame.
    % Mostly for testing the OrbitControl framework
    properties
        acc % thrust acceleration magnitude
        dir % unit vector in rsw - thrust direction
    end

    methods
        function CT = ConstantThrust(ThrustMagnitude,ThrustDirection,primary)
            %%%% Pre Initialization %%%%
            switch nargin
                case 0 % default tangential
                    ThrustMagnitude = 0.00001; % 0.001 m/s^2
                    ThrustDirection = [0,1,0].'; % s
                    primary = earth();
                case 1 % tangential
                    ThrustDirection = [0,1,0].'; % s
                    primary = earth();
                case 2 % arbitrary
                    primary = earth();
                case 3

                otherwise
                    error('Wrong number of input arguments')
            end
            %%%% Object Initialization %%%%
            % Call superclass constructor before accessing object
            CT = CT@OrbitControl(primary);
            
            %%%% Post Initialization %%%%
            % property assignment
            CT.acc = ThrustMagnitude;
            if norm(ThrustDirection) ~= 1
                ThrustDirection = ThrustDirection/norm(ThrustDirection);
            end
            CT.dir = ThrustDirection;
        end

        function [dOe] = ControlOE(CT,t,OE) %#ok<INUSD>
            % Calculates the change in Classical OE corresponding to the
            % given instantaneous orbital elements
            % handle input
            nSats = size(OE,2);
            mu = 1;
            % re = 1;
            dScale = CT.primary.Re;
            tScale = sqrt(CT.primary.Re^3/CT.primary.mu);

            normAcc = CT.acc*tScale^2/dScale;

            a = OE(1,:).';
            e = OE(2,:).';
            i = OE(3,:).';
            % ran = OE(4,:);
            w = OE(5,:).';
            M = OE(6,:).';
            f  = pi/180*me2ta(M*180/pi,e); %  conversion M to f
            p = a.*(1-e.^2);
            h = sqrt(mu*p);
            r = p./(1+e.*cos(f));
            n = sqrt(mu./a.^3);
            % eta = sqrt(1-e.^2);
            aol = f + w;
            
            % GVE
            B = zeros(nSats*6,3);

            B((0:(nSats-1))*6+1,1:2) = [2*a.^2./h.*e.*sin(f), 2*a.^2./h.*(1+e.*cos(f))]; 
            B((0:(nSats-1))*6+2,1:2) = [p./h.*sin(f), r./h.*(e+2*cos(f)+e.*cos(f).^2)];
            B((0:(nSats-1))*6+3,3) = r./h.*cos(aol);
            B((0:(nSats-1))*6+4,3) = r./h.*sin(aol)./sin(i);
            B((0:(nSats-1))*6+5,1:3) = [-p./h./e.*cos(f), r./h./e.*(2+e.*cos(f)).*sin(f),...
                -r.*sin(aol).*cos(i)./h./sin(i)];
            B((0:(nSats-1))*6+6,1:2) = [(p.*cos(f)-2*r.*e)./n./a.^2./e, -(p+r).*sin(f)./n./a.^2./e];
            % assign
            dOe = B*CT.dir*normAcc;
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
