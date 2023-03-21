classdef LatticeConstellation < Constellation
    %LatticeConstellation Defines a 3D Lattice Flower Constellation
% Based on the work in J. J. Davis, · Martín, E. Avendaño, and D. Mortari,
% “The 3-D lattice theory of Flower Constellations,” , 2013.

% To get a repeat ground track the orbital elements must be chosen
% according to the RGT condition. In addition the phasing parameters must
% fulfill:
%              (1)   nRepeats = l*gcd(nPlanes,nC3) + m*nC1
%              (2)   nDays    = m*nSatsPerAop
% where l and m are any integers.
    
    properties (SetAccess = private)
        % General Architecture
        % nSats (in superclass) - Derived
        nPlanes     % N_o
        nAops       % N_w 
        nSatsPerAop % N_so'
%         nRepeats    % N_p
%         nDays       % N_d
        
        % Phasing
        nC1     % 1-N_o
        nC2     % 1-N_w
        nC3     % 1-N_o
        
        % Orbits
        inc
        ecc
        sma
        
        % Initial Conditions
        M1
        raan1
        aop1
        
    end
    
    methods
        function LC = LatticeConstellation(Arch,Phase,Orbit,InitCon)
            %%%% Pre Initialization %%%%
            switch nargin
                case 0
                    Arch.nPlanes     = 3; % N_o
                    Arch.nAops       = 2; % N_w
                    Arch.nSatsPerAop = 1; % N'_so
                    
                    Phase.nC1    = 1;
                    Phase.nC2    = 2;
                    Phase.nC3    = 3;
                    
%                     Orbit.inc   = asind(sqrt(4/5));
                    Orbit.inc   = 50;
                    Orbit.ecc   = 0.3;
                    Orbit.sma   = CalcRgtSma(Orbit.ecc,Orbit.inc,7,1);
                    
                    InitCon.M1    = 0;
                    InitCon.raan1 = 0;
                    InitCon.aop1  = 90;
                    primary = earth();
                case 3
                    primary = earth();
                    InitCon.M1    = 0;
                    InitCon.raan1 = 0;
                    InitCon.aop1  = 0;
                case 4
                    primary = earth();
            end
            Arch.nSats = Arch.nPlanes*Arch.nAops*Arch.nSatsPerAop;
                        
            %%%% Object Initialization %%%%
            % Call superclass constructor before accessing object
            LC = LC@Constellation(Arch.nSats,primary);
            
            %%%% Post Initialization %%%%
            % property assignment
            LC.nPlanes     = Arch.nPlanes;
            LC.nAops       = Arch.nAops;
            LC.nSatsPerAop = Arch.nSatsPerAop;
            
            LC.nC1 = Phase.nC1;
            LC.nC2 = Phase.nC2;
            LC.nC3 = Phase.nC3;
            
            LC.inc = Orbit.inc;
            LC.ecc = Orbit.ecc;
            LC.sma = Orbit.sma;
            
            LC.M1    = InitCon.M1;
            LC.raan1 = InitCon.raan1;
            LC.aop1  = InitCon.aop1;
        end
        
        function oe = InitialOeOsc(LC) %[a e i O w M]
            % Returns Orbital Elements of constellation as 6xnSats matrix
            % of column vectors in the order:[a e i O w M].'
            oeM = LC.InitialOeMean();
            oe  = me2osc(oeM,LC.J2,LC.Re);
        end
        
        function X = InitialStateEci(LC)
            % Returns ECI state of constellation as 6xnSats matrix of
            % column vectors
            oe = LC.InitialOeOsc;
            oe(6,:) = me2ta(oe(6,:),oe(2,:));
            [R, V] = oe2eci(oe,LC.primary);
            X = [R; V];
        end
        
        function oeM = InitialOeMean(LC)
            % returns the orbital elements as a matrix of column vectors
            % each column represents [a,e,i,RAAN,AOP,Me].'
            X = zeros(6,LC.nSats);
            X(1,:) = LC.sma;
            X(2,:) = LC.ecc;
            X(3,:) = LC.inc;
            
            E = [LC.nPlanes,        0,              0;
                     LC.nC3, LC.nAops,              0;
                     LC.nC1,   LC.nC2, LC.nSatsPerAop];
            for iPlane = 1:LC.nPlanes
                for kAop = 1:LC.nAops
                    for jSat = 1:LC.nSatsPerAop
                        X(4:6,(iPlane-1)*LC.nSats/LC.nPlanes + ...
                            (kAop-1)*LC.nSatsPerAop + jSat) =...
                            wrapTo360(E\(360*[iPlane-1;kAop-1;jSat-1]));
                    end
                end
            end
            X(4,:) = wrapTo360(X(4,:) + LC.raan1);
            X(5,:) = wrapTo360(X(5,:) + LC.aop1);
            X(6,:) = wrapTo360(X(6,:) + LC.M1);
            oeM = X;
        end
    end
    
    methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(WC) %#ok
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Constellation Structure';
            propgroups(1).PropertyList = {'nSats','nPlanes','nAops','nSatsPerAop'};
            propgroups(2).Title = 'Constellation Phasing Parameters';
            propgroups(2).PropertyList = {'nC1','nC2','nC3'};
            propgroups(3).Title = 'Orbit Parameters';
            propgroups(3).PropertyList = {'sma','ecc','inc'};
            propgroups(4).Title = 'Primary Body Characteristics';
            propgroups(4).PropertyList = {'mu','Re','J2','wE'};
        end
    end
    
end

