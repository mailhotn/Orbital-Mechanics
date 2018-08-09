classdef LatticeConstellation < Constellation
    %LatticeConstellation Defines a 3D Lattice Flower Constellation
% Based on the work in J. J. Davis, · Martín, E. Avendaño, and D. Mortari,
% “The 3-D lattice theory of Flower Constellations,” , 2013.
% This implementation creates a constellation with a repeating earth ground
% track.
    
    properties (SetAccess = private)
        % General Architecture
        % nSats (in superclass)
        nPlanes     % N_o
        nAops       % N_w 
        nSatsPerAop % N_so'
        nRepeats    % N_p
        nDays       % N_d
        
        % Phasing
        nC1
        nC2
        nC3
        
        % Orbits
        inc
        ecc
        
        % Initial Conditions
        M1
        raan1
        aop1
        
        % Derived
        sma
    end
    
    methods
        function LC = LatticeConstellation(Arch,Phase,Orbit,Init)
            %%%% Pre Initialization %%%%
            switch nargin
                case 0
                    Arch.nSats    = 6;
                    Arch.nPlanes  = 3;
                    Arch.nAops    = 2;
                    Arch.nRepeats = 7;
                    Arch.nDays    = 1;
                    
                    Phase.nC1    = 1;
                    Phase.nC2    = 2;
                    Phase.nC3    = 3;
                    
%                     Orbit.inc   = asind(sqrt(4/5));
                    Orbit.inc   = 50;
                    Orbit.ecc   = 0.3;
                    
                    Init.M1    = 0;
                    Init.raan1 = 0;
                    Init.aop1  = 0;
                    primary = earth();
                case 3
                case 4
            end
            % Check input
            
            %%%% Object Initialization %%%%
            % Call superclass constructor before accessing object
            LC = LC@Constellation(Arch.nSats,primary);
            
            %%%% Post Initialization %%%%
            % property assignment
            LC.nPlanes     = Arch.nPlanes;
            LC.nAops       = Arch.nAops;
            LC.nRepeats    = Arch.nRepeats;
            LC.nDays       = Arch.nDays;
            LC.nSatsPerAop = LC.nSats/LC.nPlanes/LC.nAops;
            
            LC.nC1 = Phase.nC1;
            LC.nC2 = Phase.nC2;
            LC.nC3 = Phase.nC3;
            
            LC.inc = Orbit.inc;
            LC.ecc = Orbit.ecc;
            
            LC.M1    = Init.M1;
            LC.raan1 = Init.raan1;
            LC.aop1  = Init.aop1;
            
            LC.sma = CalcRgtSma(LC.ecc,LC.inc,LC.nRepeats,LC.nDays,LC.primary);
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
            [R, V] = oe2eci(oe,LC.mu);
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
            X(4,:) = X(4,:) + LC.raan1;
            X(5,:) = X(5,:) + LC.aop1;
            X(6,:) = X(6,:) + LC.M1;
            oeM = X;
        end
    end
    
    methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(WC) %#ok
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Constellation Structure';
            propgroups(1).PropertyList = {'nSats','nRepeats','nDays'};
            propgroups(2).Title = 'Constellation Phasing Parameters';
            propgroups(2).PropertyList = {'nC1','nC2','nC3'};
            propgroups(3).Title = 'Orbit Parameters';
            propgroups(3).PropertyList = {'sma','ecc','inc'};
            propgroups(4).Title = 'Primary Body Characteristics';
            propgroups(4).PropertyList = {'mu','Re','J2','wE'};
        end
    end
    
end

