classdef WalkerConstellation  < handle
    % WalkerConstellation is a class that defines an i-T/P/F Walker
    % constellation
    properties
        % constellation properties
        T   = []; % total number of satellites
        P   = []; % number of planes
        F   = []; % between-plane phasing
        inc = []; % inclination [deg]
        alt = []; % altitude [km]
        S   = []; % satellites per plane
        PU  = []; % pattern unit
        
        % environment properties
        mu  = []; % gravitational const [km^3/s^2]
        Re  = []; % primary body equatorial radius [km]
        J2  = []; % zonal harmonic
        
    end
    
    methods
        function WC = WalkerConstellation(varargin)
            switch nargin
                case 0
                    WC.T   = 15;
                    WC.P   = 3;
                    WC.inc = 55;
                    WC.alt = 1000;
                case 4 % Earth orbit
                    % constellation definition
                    WC.T   = varargin{1};
                    WC.P   = varargin{2};
                    WC.F   = varargin{3};
                    WC.inc = varargin{4};
                    WC.alt = varargin{5};
                    % primary definition
                    WC.mu  = 398600.440;
                    WC.Re  = 6378.14;
                    WC.J2  = 0.0010826265;
            end
            % derived properties
            WC.S  = WC.T/WC.P;
            WC.PU = 360/WC.T;
        end
        
        function OE = getInitElements(WC)
            % returns the orbital elements as a matrix of column vectors
            % each column represents [a,e,i,RAAN,th*]
            OE = zeros(5,WC.T);
            OE(1,:) = WC.alt + WC.Re;
            OE(3,:) = WC.inc;
            for ii = 1:WC.P
                OE(4,ii:ii+WC.S-1) = wrapTo360((ii-1)*WC.S*WC.PU);
                OE(5,ii:ii+WC.S-1) = wrapTo360((ii-1)*WC.PU*WC.F:WC.PU*WC.P
                    
        end
    end 
end