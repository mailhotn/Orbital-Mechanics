classdef WalkerConstellation < Constellation
    % WalkerConstellation is a subclass that defines an i:T/P/F Walker
    % constellation
    properties (SetAccess = private)% constellation properties
        F    % between-plane phasing
        inc  % inclination [deg]
        alt  % altitude [km]        
        S    % satellites per plane
        PU   % pattern unit
    end
    
    methods
        function WC = WalkerConstellation(T,P,F,inc,alt,primary)
            %%%% Pre Initialization %%%%
            switch nargin
                case 0 % GPS constellation in earth orbit
                    T        = 24;
                    P        = 6;
                    F        = 2;
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
            if mod(T,P) ~= 0
                error('T/P must be a whole number')
            end
            if F > P-1
                error('F must be less than P')
            end
            
            %%%% Object Initialization %%%%
            % Call superclass constructor before accessing object
            WC = WC@Constellation(T,P,primary);
            
            %%%% Post Initialization %%%%
            % property assignment
            WC.F   = F;
            WC.inc = inc;
            WC.alt = alt;
            % derived properties
            WC.S   = WC.N_sats/WC.N_planes;
            WC.PU  = 360/WC.N_sats;
        end
        
        function OE = getInitElements(WC)
            % returns the orbital elements as a matrix of column vectors
            % each column represents [a,e,i,RAAN,AOP,th*/th/Me].'
            X = zeros(6,WC.N_sats);
            X(1,:) = WC.alt + WC.R;
            X(3,:) = WC.inc;
            for ii = 1:WC.N_planes
                X(4,((ii-1)*WC.S+1):ii*WC.S) = wrapTo360((ii-1)*WC.S*WC.PU);
                X(6,((ii-1)*WC.S+1):ii*WC.S) = wrapTo360((ii-1)*WC.PU*WC.F...
                    :WC.PU*WC.N_planes:(ii-1)*WC.PU*WC.F+(WC.S-1)*WC.PU*WC.N_planes);
            end
            OE = X;
        end
        
        function X = getInitECI(WC)
            OE = WC.getInitElements;
            [R, V] = oe2eci(OE,WC.mu);
            X = [R; V];
        end
        
        function OE_m = getInitMeanElements(WC)
            
        end
    end
    
    methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(WC) %#ok
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Constellation Definitions';
            propgroups(1).PropertyList = {'N_sats','N_planes','F','inc','alt','PU','S'};
            propgroups(2).Title = 'Primary Body Characteristics';
            propgroups(2).PropertyList = {'mu','R','J2'};
        end
    end
end