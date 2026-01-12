classdef Constellation < handle  &  matlab.mixin.CustomDisplay
    % Constellation is a superclass that defines a group of satellites to be
    % propagated together.
    
    properties (SetAccess = protected)
        nSats   % Number of Satellites
        % environment
        mu % gravitational constant [km^3/s^2]
        Re % primary radius [km]
        J2 % primary J2 harmonic
        primary % More general - should use this
        third % Third Bodies
        epoch % Start epoch in Julian date
        x3AtEpoch % 6xn3 matrix of third body states at Epoch
    end
    
    methods (Access = protected)
        % Access protected as a Constellation without a type is
        % meaningless. Constructor must be called only by subclasses,
        % never independently.
        function C = Constellation(nSats, primary, third, epoch)
            switch nargin
                case 0 % single satellite constellation
                    nSats   = 1;
                    primary = Earth;
                    third  = {Moon,Sun};
                    epoch = juliandate(2000,1,1);
                case 1 % Earth orbit
                    primary = Earth;
                    third = {Moon,Sun};
                    epoch = juliandate(2000,1,1);
                case 2 % Arbitrary, no third
                    third = {};
                    epoch = juliandate(2000,1,1);
                case 3 % arbitrary bodies, default epoch
                    epoch = juliandate(2000,1,1);
                case 4 % arbitrary 
                    
                otherwise
                    error('Wrong number of input arguments')
            end
            C.nSats    = nSats;
            C.mu       = primary.mu;
            C.Re       = primary.Re;
            C.J2       = primary.J2;
            C.primary  = primary;
            C.third = third;
            C.epoch = epoch;
            C.x3AtEpoch = nan(6,length(third));
            for i3 = 1:length(third)
                [r3,v3] = planetEphemeris(epoch,primary.name,third{i3}.name);
                C.x3AtEpoch(:,i3) = [r3.';v3.'];
            end
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