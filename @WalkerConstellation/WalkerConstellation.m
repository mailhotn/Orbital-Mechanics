classdef WalkerConstellation < Constellation
    % WalkerConstellation is a subclass that defines an i:T/P/F Walker
    % constellation
    properties (SetAccess = private)% constellation properties
        phasingF      % between-plane phasing
        inc           % inclination [deg]
        alt           % altitude [km]        
        S             % satellites per plane
        PU            % pattern unit
    end
    
    methods
        function WC = WalkerConstellation(nSatsT,nPlanesP,phasingF,inc,alt,primary)
            %%%% Pre Initialization %%%%
            switch nargin
                case 0 % GPS constellation in earth orbit
                    nSatsT   = 24;
                    nPlanesP = 6;
                    phasingF = 2;
                    inc      = 55;
                    alt      = 20180;
                    primary  = earth();
                case 5 % earth orbit
                    primary  = earth();
                case 6
                    % do nothing
                otherwise
                    error('Wrong number of input arguments')
            end
            % Check input
            if mod(nSatsT,nPlanesP) ~= 0
                error(['T/P must be a whole number, T = '...
                    num2str(nSatsT) ' P = ' num2str(nPlanesP)])
            end
            if phasingF > nPlanesP-1
                error(['F must be less than P, P = ' num2str(nPlanesP) ...
                    ' F = ' num2str(phasingF)])
            end
            
            %%%% Object Initialization %%%%
            % Call superclass constructor before accessing object
            WC = WC@Constellation(nSatsT ,nPlanesP ,primary);
            
            %%%% Post Initialization %%%%
            % property assignment
            WC.phasingF = phasingF;
            WC.inc = inc;
            WC.alt = alt;
            % derived properties
            WC.S  = WC.nSats/WC.nPlanes;
            WC.PU = 360/WC.nSats;
        end
        
        function oe = InitialOeOsc(WC) %[a e i O w M]
            % Returns Orbital Elements of constealltion as 6xnSats matrix
            % of column vectors in the order:[a e i O w M].'
            oeM = WC.InitialOeMean();
            oe  = me2osc(oeM,WC.J2,WC.Re);
        end
        
        function X = InitialStateEci(WC)
            % Returns ECI state of constellation as 6xnSats matrix of
            % column vectors
            oe = WC.InitialOeOsc;
            oe(6,:) = me2ta(oe(6,:),oe(2,:));
            [R, V] = oe2eci(oe,WC.mu);
            X = [R; V];
        end
        
        function oeM = InitialOeMean(WC)
            % returns the orbital elements as a matrix of column vectors
            % each column represents [a,e,i,RAAN,AOP,Me].'
            X = zeros(6,WC.nSats);
            X(1,:) = WC.alt + WC.Re;
            X(3,:) = WC.inc;
            for ii = 1:WC.nPlanes
                X(4,((ii-1)*WC.S+1):ii*WC.S) = wrapTo360((ii-1)*WC.S*WC.PU);
                X(6,((ii-1)*WC.S+1):ii*WC.S) = wrapTo360((ii-1)*WC.PU*WC.phasingF...
                    :WC.PU*WC.nPlanes:(ii-1)*WC.PU*WC.phasingF + ...
                    (WC.S-1)*WC.PU*WC.nPlanes);
            end
            oeM = X;
        end
    end
    
    methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(WC) %#ok
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Constellation Definitions';
            propgroups(1).PropertyList = {'nSats','nPlanes','phasingF','inc',...
                                          'alt','PU','S'};
            propgroups(2).Title = 'Primary Body Characteristics';
            propgroups(2).PropertyList = {'mu','Re','J2'};
        end
    end
end