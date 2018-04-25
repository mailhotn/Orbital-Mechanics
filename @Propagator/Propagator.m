classdef Propagator < handle & matlab.mixin.Copyable & matlab.mixin.CustomDisplay
    properties(SetAccess = protected) % propagation
        reltol
        abstol
        Con     % Constellation object to propagate
    end
    
    properties(SetAccess = protected) % environment
        mu % gravitational constant [km^3/s^2]
        R  % primary radius [km]
        J2 % primary J2 harmonic
    end
    
    methods
        function P = Propagator(con, reltol, abstol, primary)
            switch nargin
                case 0 % default tolerances, earth orbit, GPS constellation
                    reltol     = 1e-8;
                    abstol     = 1e-9;
                    primary.mu = 398600.440;
                    primary.R  = 6378.137;
                    primary.J2 = 0.0010826265;
                    con = WalkerConstellation;
                case 1 % default tolerances, earth orbit
                    reltol     = 1e-8;
                    abstol     = 1e-9;
                    primary.mu = 398600.440;
                    primary.R  = 6378.137;
                    primary.J2 = 0.0010826265;
                case 2 % earth orbit
                    primary.mu = 398600.440;
                    primary.R  = 6378.137;
                    primary.J2 = 0.0010826265;
                case 3 % arbitrary orbit & tolerances
                                        
                otherwise
                    error('Wrong number of input arguments')
            end
            P.reltol = reltol;
            P.abstol = abstol;
            P.mu     = primary.mu;
            P.R      = primary.R;
            P.J2     = primary.J2;
            P.Con    = con;
        end
        
        function [Time, X] = prop_ECI_J2(P,T)
            opts = odeset('reltol',P.reltol,'abstol',P.abstol);
            [Time, X] = ode45(@P.dyn_ECI_J2,T,P.Con.getInitState,opts);
        end
    end
    
    methods(Access = protected)
        % equations of motion functions should only be called by prop
        % functions, thus they are protected
        function dX = dyn_ECI_J2(P,t,X)
            order = 6*P.Con.N_sats;
            X2 = reshape(X,[6,P.Con.N_sats]);
            Rv    = [1 0 0 0 0 0;
                     0 1 0 0 0 0;
                     0 0 1 0 0 0;
                     0 0 0 0 0 0;
                     0 0 0 0 0 0;
                     0 0 0 0 0 0]*X2;
            r_vec = sqrt(dot(Rv,Rv,1));
            r_vec = reshape([repmat(r_vec,3,1);zeros(3,P.Con.N_sats)],order,1);
            dX = repmat([1 1 1 0 0 0].',P.Con.N_sats,1) - P.mu*
        end
    end
    methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(P)
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Propagation Definitions';
            propgroups(1).PropertyList = {'reltol','abstol','Con'};
            propgroups(2).Title = 'Primary Body Characteristics';
            propgroups(2).PropertyList = {'mu','R','J2'};
        end
    end
end