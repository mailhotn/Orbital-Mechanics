classdef FlowerConstellation < Constellation
    % FlowerConstellation subclass that defines a Flower Constellation
    % nPetals-nDays-nSats-fN-fD-fH, (w-i-hP) flower constellation about the Earth.
    % The derivation is based off of D. Mortari and M. P. Wilkins, 
    % “Flower constellation set theory part I: Compatibility and phasing,” , 2008.
    % These satellites have the same repeating ground track 

    properties (SetAccess = private)% constellation properties
        % Architecture
        nPetals  % number of petals, also number of Earth-Repeats
      % nSats    % in superclass
        nDays    % number of Days to repeat
        
        % Phasing
        fN       % Phasing Parameter
        fD       % Total number of distinct orbits
        fH       % Phasing Step
        
        % Orbits
        aop      % Argument of perigee [deg]
        inc      % Inclination [deg]
        altP     % Perigee Altitude [km]
        
        % Initial Conditions
        raan0    % RAAN of first satellite [deg]
        M0       % Mean anomaly of first satellite [deg]
        
        % derived
        sma      % Semimajor Axis [km]
        raanRate % Rate of change of Mean RAAN [rad/sec]
        meanRate % Rate of change of Mean Mean Anomaly [rad/sec]
        % environment
        wE       % Earth rotation rate [rad/sec]
    end
    
    methods
        function FC = FlowerConstellation(Architecture, Phasing, Orbits, Initial)
            %%%% Pre Initialization %%%%
            switch nargin
                case 0
                    % Default Constellation
                    Architecture.nPetals = 7;
                    Architecture.nDays   = 1;
                    Architecture.nSats   = 6;
                    
                    Phasing.fD = 7;
                    Phasing.fN = 5;
                    Phasing.fH = 0;
                    
                    Orbits.inc  = 50;
                    Orbits.aop  = 0;
                    Orbits.altP = 1.6751e+03;
                    
                    Initial.raan0 = 0;
                    Initial.M0 = 0;
                    
                    primary = earth();
                case 3
                    Initial.raan0 = 0;
                    Initial.M0 = 0;
                    primary = earth();
                case 4
                    primary = earth();         
                
                otherwise
                    error('Wrong number of input arguments')
            end
            % Check input
            if (gcd(Architecture.nPetals,Architecture.nDays) ~= 1) && ...
                    (Architecture.nPetals ~= Architecture.nDays)
                Ct = gcd(Architecture.nPetals,Architecture.nDays);
                Architecture.nPetals = Architecture.nPetals/Ct;
                Architecture.nDays = Architecture.nDays/Ct;
            end
            if (Phasing.fN > Phasing.fD)
                Phasing.fN = mod(Phasing.fN,Phasing.fD);
            end
            if (gcd(Phasing.fN,Phasing.fD) ~= 1) && ...
                    (Phasing.fN ~= Phasing.fD)
                error('fN & fD are not relatively prime');
            end
            if (Phasing.fH > (Architecture.nDays - 1))
                error('fH is too large');
            end
            
            %%%% Object Initialization %%%%
            % Call superclass constructor before accessing object
            FC = FC@Constellation(Architecture.nSats,primary);
            
            %%%% Post Initialization %%%%
            % property assignment
            FC.wE = primary.we*pi/180;
            
            FC.nPetals = Architecture.nPetals;
            FC.nDays   = Architecture.nDays;
            
            FC.fN = Phasing.fN;
            FC.fD = Phasing.fD;
            FC.fH = Phasing.fH;
            
            FC.aop  = Orbits.aop;
            FC.inc  = Orbits.inc;
            FC.altP = Orbits.altP;
            
            FC.raan0 = Initial.raan0;
            FC.M0 = Initial.M0;
            
            FC.CalcRgtOrbit();
        end
                
        function oe = InitialOeOsc(FC) %[a e i O w M]
            % Returns Orbital Elements of constellation as 6xnSats matrix
            % of column vectors in the order:[a e i O w M].'
            oeM = FC.InitialOeMean();
            oe  = me2osc(oeM,FC.J2,FC.Re);
        end
        
        function X = InitialStateEci(FC)
            % Returns ECI state of constellation as 6xnSats matrix of
            % column vectors
            oe = FC.InitialOeOsc;
            oe(6,:) = me2ta(oe(6,:),oe(2,:));
            [R, V] = oe2eci(oe,FC.mu);
            X = [R; V];
        end
        
        function oeM = InitialOeMean(FC)
            % returns the orbital elements as a matrix of column vectors
            % each column represents [a,e,i,RAAN,AOP,Me].'
            n   = (FC.sma^3/FC.mu)^-(1/2);
            X = zeros(6,FC.nSats);
            X(1,:) = FC.sma;
            X(2,:) = 1 - (FC.Re + FC.altP)/FC.sma;
            X(3,:) = FC.inc;
            X(4,1) = FC.raan0;
            X(5,:) = FC.aop;
            X(6,1) = FC.M0;
            % Assign RAAN
            for iSat = 1:(FC.nSats-1)
                X(4,iSat+1) = X(4,iSat) + 360*FC.fN/FC.fD;
            end
            % Assign Mean Anomaly
            for iSat = 1:(FC.nSats-1)
                X(6,iSat+1) =  X(6,iSat) - ...
                    360*FC.fN/FC.fD*(n+FC.meanRate)/(FC.wE - FC.raanRate) - ...
                    360*FC.fH/FC.nDays;
            end
            oeM = X;
        end
    end
    
    methods(Access = protected)
        % Internal  use only, usesd to assign properties in object
        % creation.
        function  CalcRgtOrbit(FC)
            % N-R setup
            tol = 1e-12;
            ratio = inf;
            
            incRad = FC.inc*pi/180;
            tau = FC.nDays/FC.nPetals;
            % initial guess
            a = ((tau/FC.wE)^2*FC.mu)^(1/3);
            % Newton Raphson Loop
            while abs(ratio) > tol
                p   = 2*(FC.Re + FC.altP) - (FC.Re + FC.altP)^2/a;
                e   = 1 - (FC.Re + FC.altP)/a;
                xi  = 3*FC.Re^2*FC.J2/(4*p^2);
                chi = 1 + xi*(4 + 2*sqrt(1-e^2) - (5 + 3*sqrt(1-e^2)*(sin(incRad))^2));
                n   = (a^3/FC.mu)^-(1/2);
                
                f = FC.mu^(-1/2)*a^(3/2) - tau/FC.wE*chi/(1 + 2/FC.wE*xi*n*cos(incRad));
                df = (3*a^(1/2))/(2*FC.mu^(1/2)) - (tau*((3*FC.J2*FC.Re^2*(4*a + 2*sign(a)*(-(FC.Re + FC.altP)*(FC.Re - 2*a + FC.altP))^(1/2) - 5*a*sin(incRad)^2 - 3*sign(a)*sin(incRad)^2*(-(FC.Re + FC.altP)*(FC.Re - 2*a + FC.altP))^(1/2)))/(2*(FC.Re + FC.altP)*(FC.Re - 2*a + FC.altP)^3) - (3*FC.J2*FC.Re^2*sign(a)*(3*sin(incRad)^2 - 2)*(FC.Re - a + FC.altP))/(4*(FC.Re + FC.altP)*(-(FC.Re + FC.altP)*(FC.Re - 2*a + FC.altP))^(1/2)*(FC.Re - 2*a + FC.altP)^2)))/(FC.wE*((3*FC.J2*FC.Re^2*a^2*cos(incRad)*(FC.mu/a^3)^(1/2))/(2*FC.wE*(FC.Re + FC.altP)^2*(FC.Re - 2*a + FC.altP)^2) + 1)) + (3*FC.J2*FC.Re^2*FC.mu*tau*cos(incRad)*((3*FC.J2*FC.Re^2*a*(4*a + 2*sign(a)*(-(FC.Re + FC.altP)*(FC.Re - 2*a + FC.altP))^(1/2) - 5*a*sin(incRad)^2 - 3*sign(a)*sin(incRad)^2*(-(FC.Re + FC.altP)*(FC.Re - 2*a + FC.altP))^(1/2)))/(4*(FC.Re + FC.altP)^2*(FC.Re - 2*a + FC.altP)^2) + 1)*(FC.Re + 6*a + FC.altP))/(4*a^2*FC.wE^2*(FC.Re + FC.altP)^2*(FC.mu/a^3)^(1/2)*((3*FC.J2*FC.Re^2*a^2*cos(incRad)*(FC.mu/a^3)^(1/2))/(2*FC.wE*(FC.altP + FC.Re)^2*(FC.altP - 2*a + FC.Re)^2) + 1)^2*(FC.Re - 2*a + FC.altP)^3);
                ratio = f/df;
                a = a - ratio;
            end
            % Re-evaluate
            p   = 2*(FC.Re + FC.altP) - (FC.Re + FC.altP)^2/a;
            e   = 1 - (FC.Re + FC.altP)/a;
            xi  = 3*FC.Re^2*FC.J2/(4*p^2);
            n   = (a^3/FC.mu)^-(1/2);
            % Assign Values
            FC.sma = a;
            FC.raanRate = -2*xi*n*cos(incRad);
            FC.meanRate = -xi*n*sqrt(1-e^2)*(3*sin(incRad)^2 - 2);
        end
    end
    
    methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(WC) %#ok
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Constellation Structure';
            propgroups(1).PropertyList = {'nPetals','nSats','nDays'};
            propgroups(2).Title = 'Constellation Phasing Parameters';
            propgroups(2).PropertyList = {'fD','fN','fH'};
            propgroups(3).Title = 'Orbit Parameters';
            propgroups(3).PropertyList = {'sma','altP','inc','aop'};
            propgroups(4).Title = 'Primary Body Characteristics';
            propgroups(4).PropertyList = {'mu','Re','J2','wE'};
        end
    end
end