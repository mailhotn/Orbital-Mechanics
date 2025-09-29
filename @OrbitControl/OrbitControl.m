classdef OrbitControl < handle & matlab.mixin.CustomDisplay
    % OrbitControl is a superclass that defines a controller for use with
    % the Propagator class.

    properties (SetAccess = protected)
        primary % environment
    end

    methods (Access = protected)
        % Access protected as a Controller without a type is
        % meaningless. Constructor must be called only by subclasses,
        % never independently.
        function C = OrbitControl(primary)
            switch nargin
                case 0
                    primary = earth();
                case 1 % Arbitrary

                otherwise
                    error('Wrong number of input arguments')
            end
            C.primary  = primary;
        end
    end
end