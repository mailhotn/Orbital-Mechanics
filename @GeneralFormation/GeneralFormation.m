classdef GeneralFormation < Constellation
    %GeneralFormation Defines a chief and any number of deputy satellites.
% The formation is defined through the chief's orbital elements and the
% deputies' differential orbital elements

    properties (SetAccess = private)
        % General Architecture
        % nSats (in superclass) - Derived
        
        % Chief Elements
        oeChief
        
        % Deputy Differential Elements
        dOeDeputy
        
        % Brouwer Flag (For DSS Homework)
        useBrouwer
        
    end
    
    methods
        function GF = GeneralFormation(oeChief,dOeDeputy,brouwer,primary)
            %%%% Pre Initialization %%%%
            switch nargin
                case 0
                    oeChief = [7100, 0.05, 90, 0, 30, 0].';
                    dOeDeputy = [0,0,0,0,0,1].';
                    brouwer = true;
                    primary = earth();
                case 2
                    brouwer = true;
                    primary = earth();
                case 3
                    primary = earth();
                case 4
                    
            end
            nSats = 1 + size(dOeDeputy,2);
                        
            %%%% Object Initialization %%%%
            % Call superclass constructor before accessing object
            GF = GF@Constellation(nSats,primary);
            
            %%%% Post Initialization %%%%
            % property assignment
            GF.oeChief   = oeChief;
            GF.dOeDeputy = dOeDeputy;

            GF.useBrouwer = brouwer;
        end
        
        function oe = InitialOeOsc(GF) %[a e i O w M]
            % Returns Orbital Elements of constellation as 6xnSats matrix
            % of column vectors in the order:[a e i O w M].'
            oeM = GF.InitialOeMean();
            if GF.useBrouwer
                oe = zeros(6,GF.nSats);
                for iSat = 1:GF.nSats
                    oe(:,iSat)  = me2oscNum(oeM(:,iSat));
                end
            else
                oe = oeM;
            end
        end
        
        function X = InitialStateEci(GF)
            % Returns ECI state of constellation as 6xnSats matrix of
            % column vectors
            oe = GF.InitialOeOsc;
            oe(6,:) = me2ta(oe(6,:),oe(2,:));
            [R, V] = oe2eci(oe,GF.primary);
            X = [R; V];
        end
        
        function oeM = InitialOeMean(GF)
            % returns the orbital elements as a matrix of column vectors
            % each column represents [a,e,i,RAAN,AOP,Me].'
            X = zeros(6,GF.nSats);
            X(:,1) = GF.oeChief;
            X(:,2:end) = GF.oeChief + GF.dOeDeputy;
            
            oeM = X;
        end

        function oeM = InitialOeMeanShort(GF) %???
            % returns the orbital elements as a matrix of column vectors
            % each column represents [a,e,i,RAAN,AOP,Me].'
            X = zeros(6,GF.nSats);
            X(:,1) = GF.oeChief;
            X(:,2:end) = GF.oeChief + GF.dOeDeputy;
            
            oeM = X;
        end
    end
    
    methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(WC) %#ok
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Chief Orbital Elements';
            propgroups(1).PropertyList = {'oeChief'};
            propgroups(2).Title = 'Deputy Differential Elements';
            propgroups(2).PropertyList = {'dOeDeputy'};
            propgroups(3).Title = 'Primary Body Characteristics';
            propgroups(3).PropertyList = {'mu','Re','J2','wE'};
        end
    end  
end