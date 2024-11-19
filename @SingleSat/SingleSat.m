classdef SingleSat < Constellation
    % SingleSat defines a "constellation" of a single satellite with given
    % osculating orbital elements. This is useful for using the constellation
    % propagator framework when working with a single satellite.
    properties
        sma  % semi-major axis [km]
        e    % eccentricity
        inc  % inclination [deg]
        raan % right ascension of ascending node [deg]
        aop  % argument of perigee [deg]
        me   % Mean anomaly [deg]
    end
    
    methods
        function Sat = SingleSat(OE,primary)
            %%%% Pre Initialization %%%%
            switch nargin
                case 0 % default - Hubble
                    OE = [6917.5, 0.000287, 28.47,...
                          176.23, 82.61, 319.41].';
                    primary = earth();
                case 1 % earth orbit
                    primary = earth();
                case 2
                    
                otherwise
                    error('Wrong number of input arguments')
            end
            %%%% Object Initialization %%%%
            % Call superclass constructor before accessing object
            Sat = Sat@Constellation(1,primary);
            
            %%%% Post Initialization %%%%
            % property assignment
            Sat.sma  = OE(1);
            Sat.e    = OE(2);
            Sat.inc  = OE(3);
            Sat.raan = OE(4);
            Sat.aop  = OE(5);
            Sat.me   = OE(6);
        end
        
        function oeM = InitialOeMean(Sat)
            oeOsc = Sat.InitialOeOsc;
            oeM = osc2me(oeOsc,Sat.J2,Sat.Re);
        end
        
        function oeOsc = InitialOeOsc(Sat)            
            oeOsc = [Sat.sma,Sat.e,Sat.inc,Sat.raan,Sat.aop,Sat.me].';
        end
        
        function X = InitialStateEci(Sat)
            oeOsc = Sat.InitialOeOsc;
            oeOsc(6) = me2ta(oeOsc(6),oeOsc(2));
            [R, V] = oe2eci(oeOsc,Sat.primary);
            X = [R; V];
        end
        
        function oeM = InitialOeMeanShort(Sat)
            oeOsc = Sat.InitialOeOsc;
            oeM = osc2meSP(oeOsc,Sat.primary);
        end
    end
    
     methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(Sat) %#ok
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Oscullating Orbital Elements';
            propgroups(1).PropertyList = {'sma','e','inc','raan','aop','me'};
            propgroups(2).Title = 'Primary Body Characteristics';
            propgroups(2).PropertyList = {'mu','R','J2'};
        end
    end
end