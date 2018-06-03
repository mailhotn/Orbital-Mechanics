classdef FlowerConstellation < Constellation
    % FlowerConstellation is a subclass that defines a 
    % nPetals-nDays-nSats-w-i-hP flower constellation.
    % The derivation is based off of "The Flower Constellations" - Daniele
    % Mortari, Matthew P. Wilkins, and Christian Bruccoleri, 2004
    % This constellation has a repeating ground track and is critically
    % inclined
    %
    % ~~~~~~~~~~ Parameters ~~~~~~~~~~
    % nPetals - Number of orbits until repeat, number of "petals"
    % nDays - Number of days to Repeat
    % 
    properties (SetAccess = private)% constellation properties

    end
    
    methods
        function FC = FlowerConstellation()
            %%%% Pre Initialization %%%%
            switch nargin
                
                otherwise
                    error('Wrong number of input arguments')
            end
            % Check input

            %%%% Object Initialization %%%%
            % Call superclass constructor before accessing object
            FC = FC@Constellation();

        end
        
        function OE = getInitElements(FC) %[a e i O w M]

        end
        
        function X = getInitECI(FC)

        end
        
        function OE_m = getInitMeanElements(FC)

        end
    end
    
    methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(FC) %#ok
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Constellation Definitions';
            propgroups(1).PropertyList = {'N_sats','N_planes','F','inc','alt','PU','S'};
            propgroups(2).Title = 'Primary Body Characteristics';
            propgroups(2).PropertyList = {'mu','Re','J2'};
        end
    end
end