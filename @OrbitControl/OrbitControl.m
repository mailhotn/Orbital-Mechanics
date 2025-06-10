classdef OrbitControl < handle & matlab.mixin.CustomDisplay
    % OrbitControl is a superclass that defines a controller for use with
    % the Propagator class.

        properties (SetAccess = protected)
        Con     % Constellation
        primary % environment
        end

        methods (Access = protected)
        % Access protected as a Controller without a type is
        % meaningless. Constructor must be called only by subclasses,
        % never independently.
        function C = OrbitControl(Constellation, primary)
            switch nargin
                case 0 % single satellite constellation
                    Constellation = SingleSat;
                    primary = earth();
                case 1 % Earth orbit
                    primary = earth();
                case 2 % Arbitrary
                    
                otherwise
                    error('Wrong number of input arguments')
            end
            C.Con      = Constellation;
            C.primary  = primary;
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