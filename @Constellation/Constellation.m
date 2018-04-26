classdef Constellation < handle  &  matlab.mixin.CustomDisplay
    % Constellation is a superclass that defines a group of satellites to be
    % propagated together.
    
    properties (SetAccess = protected)
        N_sats   % Number of Satellites
        N_planes % Number of planes
    end
    
    properties(SetAccess = protected) % environment
        mu % gravitational constant [km^3/s^2]
        R  % primary radius [km]
        J2 % primary J2 harmonic
    end
    
    methods (Access = protected)
        % Access protected as a Constellation without a type is
        % meaningless. Constructor must be called obnly by subclasses,
        % never independantly.
        function C = Constellation(N_sats, N_planes, primary)
            switch nargin
                case 0 % single satellite constellation
                    N_sats   = 1;
                    N_planes = 1;
                    primary = earth();
                case 2 % Earth orbit
                    primary = earth();
                case 3 % Arbitrary
                    
                otherwise
                    error('Wrong number of input arguments')
            end
            C.N_sats   = N_sats;
            C.N_planes = N_planes;
            C.mu       = primary.mu;
            C.R        = primary.R;
            C.J2       = primary.J2;
        end
    end
%     ~~~~~~~~~~~~~ Moved To Subclasses ~~~~~~~~~~~~~
%     methods(Access = protected)
%         function propgroups = getPropertyGroups(C)
%             propgroups = matlab.mixin.util.PropertyGroup;
%             propgroups(1).Title = 'Constellation Definitions';
%             propgroups(1).PropertyList = {'N_sats','N_planes'};
%             propgroups(2).Title = 'Primary Body Characteristics';
%             propgroups(2).PropertyList = {'mu','R','J2'};
%         end
%     end
end