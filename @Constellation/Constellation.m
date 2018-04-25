classdef Constellation < handle & matlab.mixin.Copyable
    % Constellation is a class that defines a group of satellites to be
    % propagated together.
    properties (SetAccess = protected)
        N_sats   % Number of Satellites
        N_planes % Number of planes
    end
    
    methods
        function C = Constellation(N_sats, N_planes)
            switch nargin
                case 0
                    N_sats   = 24;
                    N_planes = 6;
                case 2
                    
            end
            C.N_sats   = N_sats;
            C.N_planes = N_planes;
        end
    end
end