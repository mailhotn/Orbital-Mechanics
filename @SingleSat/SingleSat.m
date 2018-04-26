classdef SingleSat < Constellation
    % SingleSat defines a "constellation" of a single satellite with given
    % orbital elements. This is useful for using the constellation
    % propagator framework when working with a single satellite.
    properties
        sma  % semi-major axis [km
        e    % eccentricity
        inc  % inclination [deg]
        raan % right ascension of ascending node [deg]
        aop  % argument of perigee [deg]
        th   % true anomaly [deg]
    end
    
    methods
        function Sat = SingleSat(OE,primary)
            %%%% Pre Initialization %%%%
            switch nargin
                case 0 % default - Hubble
                    OE = [6917.5, 0.0002455, 28.4690,...
                          103.9640, 294.6441, 0].';
                    primary = earth();
                case 1 % earth orbit
                    primary = earth();
                case 2
                    
                otherwise
                    error('Wrong numnber of input arguments')
            end
            %%%% Object Initialization %%%%
            % Call superclass constructor before accessing object
            Sat = Sat@Constellation(1,1,primary);
            
            %%%% Post Initialization %%%%
            % property assignment
            Sat.sma  = OE(1);
            Sat.e    = OE(2);
            Sat.inc  = OE(3);
            Sat.raan = OE(4);
            Sat.aop  = OE(5);
            Sat.th   = OE(6);
        end
        
        function OE = getInitElements(Sat)
            OE = [Sat.sma,Sat.e,Sat.inc,Sat.raan,Sat.aop,Sat.th].';
        end
        
        function X = getInitECI(Sat)
            OE = [Sat.sma,Sat.e,Sat.inc,Sat.raan,Sat.aop,Sat.th].';
            [R, V] = oe2eci(OE,Sat.mu);
            X = [R; V];
        end
    end
    
     methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(Sat) %#ok
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Satellite Orbital Elements';
            propgroups(1).PropertyList = {'sma','e','inc','raan','aop','th'};
            propgroups(2).Title = 'Primary Body Characteristics';
            propgroups(2).PropertyList = {'mu','R','J2'};
        end
    end
end