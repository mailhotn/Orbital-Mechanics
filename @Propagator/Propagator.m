classdef Propagator < handle &  matlab.mixin.CustomDisplay
    % Propagator is a class that defines an orbit propagator for use with
    % a constellation object.
    properties(SetAccess = protected) % propagation
        reltol
        abstol
        Con     % Constellation object to propagate
    end
    
    methods
        function P = Propagator(con, reltol, abstol)
            switch nargin
                case 0 % default tolerances, GPS constellation
                    reltol     = 1e-8;
                    abstol     = 1e-9;
                    con = WalkerConstellation;
                case 1 % default tolerances 
                    reltol     = 1e-8;
                    abstol     = 1e-9;
                case 3 % arbitrary orbit & tolerances
                                        
                otherwise
                    error('Wrong number of input arguments')
            end
            P.reltol = reltol;
            P.abstol = abstol;
            P.Con    = con;
        end
        
        function [Time, X] = prop_ECI_TB(P,T)
            % Propagate for time T in ECI frame with no perturbations
            % Complexity ~O(T)
            opts = odeset('reltol',P.reltol,'abstol',P.abstol);
            IC = reshape(P.Con.getInitECI,[6*P.Con.N_sats,1]);
            [Time, X] = ode45(@P.dyn_ECI_TB,T,IC,opts);
        end
        
        function [Time, X] = prop_ECI_TB_for(P,T)
            % Propagate for time T in ECI frame with no perturbations
            % Much slower for large constellations. Slightly faster for
            % single satellites. Complexity ~O(N_sats*T)
            opts = odeset('reltol',P.reltol,'abstol',P.abstol);
            IC = P.Con.getInitECI;
            X = [];
            for ii = 1:P.Con.N_sats
                 [Time, tX] = ode45(@P.dyn_ECI_TB_for,T,IC(:,ii),opts);
                 X = [X,tX];%#ok
            end
        end
        
        function [Time, X] = prop_ECI_TB_for2(P,T)
            % Propagate for time T in ECI frame with no perturbations
            % Slower than vectorized, faster than regular for style
            % Much easier to write than vectorized
            % Complexity ~O(T)
            opts = odeset('reltol',P.reltol,'abstol',P.abstol);
            IC = reshape(P.Con.getInitECI,[6*P.Con.N_sats,1]);
            [Time, X] = ode45(@P.dyn_ECI_TB_for2,T,IC,opts);
        end
        
        function [Time, X] = prop_ECI_J2(P,T)
            % Propagate for time T in ECI frame with J2 perturbation
            % Complexity ~O(T)
            opts = odeset('reltol',P.reltol,'abstol',P.abstol);
            IC = reshape(P.Con.getInitECI,[6*P.Con.N_sats,1]);
            [Time, X] = ode45(@P.dyn_ECI_J2,T,IC,opts);
        end
        
        function [Time, X] = prop_OE_Mean(P,T)
            
        end
    end
    
    methods(Access = protected) % Equations of Motion
        % equations of motion functions should only be called by prop
        % functions, thus they are protected
        
        function dX = dyn_ECI_TB(P,t,X) %#ok
            % move stuff around
            order = 6*P.Con.N_sats;
            X2 = reshape(X,[6,P.Con.N_sats]);
            Rv = [eye(3),zeros(3);
                  zeros(3,6)]*X2;
            % get vector of R magnitudes
            r_vec = sqrt(dot(Rv,Rv,1));
            r_vec = reshape(repmat(r_vec,6,1),order,1);
            % move more stuff around
            R2 = repmat([1 1 1 0 0 0].',P.Con.N_sats,1).*X;
            V2 = repmat([0 0 0 1 1 1].',P.Con.N_sats,1).*X;
            % equation of motion
            dX = circshift(V2,-3) - circshift(P.Con.mu*R2./r_vec.^3,3);
        end
        
        function dX = dyn_ECI_TB_for(P,t,X) %#ok
            r = norm(X(1:3));
            dX = [X(4:6);-P.Con.mu*X(1:3)/r^3];
        end
        
        function dX = dyn_ECI_TB_for2(P,t,X) %#ok
            order = 6*P.Con.N_sats;
            dX = zeros(order,1);
            for ii = 1:P.Con.N_sats
                ind = (ii-1)*6+1:(ii-1)*6+6;
                R = X(ind(1:3));
                V = X(ind(4:6));
                r = norm(R);
                dX(ind) = [V;-P.Con.mu*R/r^3];
            end
        end
        
        function dX = dyn_ECI_J2(P,t,X) %#ok
            % move stuff around
            order = 6*P.Con.N_sats;
            X2 = reshape(X,[6,P.Con.N_sats]);
            Rv = [eye(3),zeros(3);
                zeros(3,6)]*X2;
            % get vector of R magnitudes
            r = sqrt(dot(Rv,Rv,1));
            r = reshape(repmat(r,6,1),order,1);
            % move more stuff around
            R2 = repmat([1 1 1 0 0 0].',P.Con.N_sats,1).*X;
            V2 = repmat([0 0 0 1 1 1].',P.Con.N_sats,1).*X;
            Z  = repmat([0 0 1 0 0 0].',P.Con.N_sats,1).*R2;
            Z2 = [zeros(3,2),ones(3,1),zeros(3);zeros(3,6)]*X2;
            Z2 = reshape(Z2,[order,1]);
            % equation of motion
            f_J2 = -P.Con.mu*P.Con.J2*P.Con.R^2./r.^4.*...
                   (3*Z./r + (-7.5*(Z2./r).^2 + 1.5).*R2./r);
            dX = circshift(V2,-3) + circshift(-P.Con.mu*R2./r.^3 + f_J2,3);
        end
    end
    
    methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(P) %#ok
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Propagation Definitions';
            propgroups(1).PropertyList = {'reltol','abstol','Con'};
        end
    end
end