classdef WalkerConstellation < Constellation
    % WalkerConstellation is a class that defines an i:T/P/F Walker
    % constellation
    properties (SetAccess = private)
        % constellation properties
        F    % between-plane phasing
        inc  % inclination [deg]
        sma  % semi-major axis [km]        
        S    % satellites per plane
        PU   % pattern unit
    end
    
    methods
        function WC = WalkerConstellation(T,P,F,inc,sma)
            %%%% Pre Initialization %%%%
            switch nargin
                case 0
                    % constellation definition
                    T        = 24;
                    P        = 6;
                    F        = 2;
                    inc      = 55;
                    sma      = 20180 + 6378.137;
                case 5
                    % Check input
                    if mod(T,P) ~= 0
                        error('Invalid T/P configuration')
                    end
                    if F > P-1
                        error('F must be less than P')
                    end
                otherwise
                    error('Wrong number of input arguments')
            end
            
            %%%% Object Initialization %%%%
            % Call superclass constructor before accessing object
            WC = WC@Constellation(T,P);
            
            %%%% Post Initialization %%%%
            % property assignment
            WC.F   = F;
            WC.inc = inc;
            WC.sma = sma;
            % derived properties
            WC.S   = WC.N_sats/WC.N_planes;
            WC.PU  = 360/WC.N_sats;
        end
        
        function OE = getInitElements(WC)
            % returns the orbital elements as a matrix of column vectors
            % each column represents [a,e,i,RAAN,AOP,th*/th/Me].'
            X = zeros(6,WC.N_sats);
            X(1,:) = WC.sma;
            X(3,:) = WC.inc;
            for ii = 1:WC.N_planes
                X(4,((ii-1)*WC.S+1):ii*WC.S) = wrapTo360((ii-1)*WC.S*WC.PU);
                X(6,((ii-1)*WC.S+1):ii*WC.S) = wrapTo360((ii-1)*WC.PU*WC.F...
                    :WC.PU*WC.N_planes:(ii-1)*WC.PU*WC.F+(WC.S-1)*WC.PU*WC.N_planes);
            end
            OE = X;
        end
    end
end