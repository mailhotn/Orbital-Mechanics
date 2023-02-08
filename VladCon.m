classdef VladCon < Constellation
    % VladCon defines a constellation for checking Vlad's homework
    properties
        sma       % semi-major axis [km]
        ecc       % eccentricity
        incVec    % inclination range 1xnSats [deg]
        ran       % right ascension of ascending node [deg]
        aop       % argument of perigee [deg]
        man       % true anomaly [deg]
    end
    
    methods
        function Con = VladCon(OE,inc,primary)
            %%%% Pre Initialization %%%%
            switch nargin
                case 0 % default - Hubble
                    OE = [6917.5, 0.000287, 28.47,...
                          176.23, 82.61, 319.41].';
                      inc = 28.47;
                      nSats = 1;
                    primary = earth();
                case 2 % earth orbit
                    nSats = length(inc);
                    primary = earth();
                case 3
                    nSats = length(inc);
                otherwise
                    error('Wrong number of input arguments')
            end
            %%%% Object Initialization %%%%
            % Call superclass constructor before accessing object
            Con = Con@Constellation(nSats,primary);
            
            %%%% Post Initialization %%%%
            % property assignment
            Con.sma  = OE(1);
            Con.ecc  = OE(2);
            Con.ran  = OE(4);
            Con.aop  = OE(5);
            Con.man  = ta2me(OE(6),OE(2));
            
            Con.incVec = inc;
        end
        
        function oeM = InitialOeMean(Con)
            oeOsc = Con.InitialOeOsc;
            oeM = osc2me(oeOsc,Con.J2,Con.Re);
        end
        
        function oeOsc = InitialOeOsc(Con)            
            oeOsc = repmat([Con.sma,Con.ecc,nan,Con.ran,Con.aop,Con.man].',1,Con.nSats);
            oeOsc(3,:) = Con.incVec;
        end
        
        function X = InitialStateEci(Con)
            oeOsc = Con.InitialOeOsc;
            [R, V] = oe2eci(oeOsc,Con.primary);
            X = [R; V];
        end
        
        function oeM = InitialOeMeanShort(Con)
            oeOsc = Con.InitialOeOsc;
            oeM = osc2meSP(oeOsc,Con.primary);
        end
    end
    
     methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(Sat) %#ok
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Satellite Orbital Elements';
            propgroups(1).PropertyList = {'sma','ecc','inc','ran','aop','man'};
            propgroups(2).Title = 'Primary Body Characteristics';
            propgroups(2).PropertyList = {'mu','R','J2'};
        end
    end
end