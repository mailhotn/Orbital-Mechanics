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
            % This should be replaced with Analytical Solution
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
            opts = odeset('reltol',P.reltol,'abstol',P.abstol);
            IC = reshape(P.Con.getInitMeanElements,[6*P.Con.N_sats,1]);
            [Time, X] = ode45(@P.dyn_OE_Mean,T,IC,opts);
        end
        
        function [Time, X] = prop_OE_Mean_lin(P,T)
            IC = reshape(P.Con.getInitMeanElements,[6*P.Con.N_sats,1]);
            dT = T(2) - T(1);
            % move stuff around
            order = 6*P.Con.N_sats;
            X2 = reshape(IC,[6,P.Con.N_sats]);
            aV = X2(1,:);
            eV = X2(2,:);
            iV = X2(3,:);
            a = reshape(repmat(aV,6,1),order,1);
            e = reshape(repmat(eV,6,1),order,1);
            i = reshape(repmat(iV,6,1),order,1);
            % Derived values
            p = a.*(1-e.^2);
            n = sqrt(P.Con.mu./a.^3);
            eta = sqrt(1-e.^2);
            % Eq of motion
            dO = 180/pi*repmat([0,0,0,1,0,0].',P.Con.N_sats,1).*...
                (-3/2*P.Con.J2.*(P.Con.Re./p).^2.*n.*cosd(i));
            dw = 180/pi*repmat([0,0,0,0,1,0].',P.Con.N_sats,1).*...
                (3/4*P.Con.J2.*(P.Con.Re./p).^2.*n.*(5*cosd(i).^2-1));
            dM = 180/pi*repmat([0,0,0,0,0,1].',P.Con.N_sats,1).*...
                (3/4*P.Con.J2.*(P.Con.Re./p).^2.*n.*eta.*(3*cosd(i).^2-1) + n);
            dX = dO + dw + dM;
            % Propagate
            X = zeros(order,length(T));
            X(:,1) = IC;
            X(:,2:end) = X(:,1) + dX*T(2:end);
            X = X.';
            Time = T;
        end
        
        function [Time, X] = prop_OE_Osc(P,T)
            opts = odeset('reltol',P.reltol,'abstol',P.abstol);
            OE = P.Con.getInitElements;
            OE(3:end,:) = OE(3:end,:)*pi/180;
            IC = reshape(OE,[6*P.Con.N_sats,1]);
            [Time, X] = ode45(@P.dyn_OE_Osc,T,IC,opts);
            inddeg = reshape((3:6).'+(0:P.Con.N_sats-1)*6,4*P.Con.N_sats,1);
            X(:,inddeg) = wrapTo360(180/pi*X(:,inddeg));
        end
    end
    
    methods(Access = protected) % Equations of Motion
        % equations of motion functions should only be called by prop
        % functions, thus they are protected
        
        function dX = dyn_ECI_TB(P,t,X)  %#ok<INUSL>
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
        
        function dX = dyn_ECI_TB_for(P,t,X) %#ok<INUSL>
            r = norm(X(1:3));
            dX = [X(4:6);-P.Con.mu*X(1:3)/r^3];
        end
        
        function dX = dyn_ECI_TB_for2(P,t,X) %#ok<INUSL>
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
        
        function dX = dyn_ECI_J2(P,t,X) %#ok<INUSL>
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
            f_J2 = -P.Con.mu*P.Con.J2*P.Con.Re^2./r.^4.*...
                   (3*Z./r + (-7.5*(Z2./r).^2 + 1.5).*R2./r);
            dX = circshift(V2,-3) + circshift(-P.Con.mu*R2./r.^3 + f_J2,3);
        end
        
        function dX = dyn_OE_Mean(P,t,X) %#ok<INUSL>
            % move stuff around
            order = 6*P.Con.N_sats;
            X2 = reshape(X,[6,P.Con.N_sats]);
            aV = X2(1,:);
            eV = X2(2,:);
            iV = X2(3,:);
            a = reshape(repmat(aV,6,1),order,1);
            e = reshape(repmat(eV,6,1),order,1);
            i = reshape(repmat(iV,6,1),order,1);
            % Derived values
            p = a.*(1-e.^2);
            n = sqrt(P.Con.mu./a.^3);
            eta = sqrt(1-e.^2);
            % Eq of motion
            dO = 180/pi*repmat([0,0,0,1,0,0].',P.Con.N_sats,1).*...
                (-3/2*P.Con.J2.*(P.Con.Re./p).^2.*n.*cosd(i));
            dw = 180/pi*repmat([0,0,0,0,1,0].',P.Con.N_sats,1).*...
                (3/4*P.Con.J2.*(P.Con.Re./p).^2.*n.*(5*cosd(i).^2-1));
            dM = 180/pi*repmat([0,0,0,0,0,1].',P.Con.N_sats,1).*...
                (3/4*P.Con.J2.*(P.Con.Re./p).^2.*n.*eta.*(3*cosd(i).^2-1) + n);
            dX = dO + dw + dM;
        end
        
        function dX = dyn_OE_Osc(P,t,X)  %#ok<INUSL>
            OE  = reshape(X,6,P.Con.N_sats);
            % Element vectors (angles already in radians)
            a   = OE(1,:);
            e   = OE(2,:);
            inc = OE(3,:);
            w   = OE(5,:);
            Me  = OE(6,:);
            th  = pi/180*me2ta(Me*180/pi,e);
            % Secondary definitions
            p = a.*(1-e.^2);
            h = sqrt(P.Con.mu*p);
            r = p./(1+e.*cos(th));
            MJ2 = -3*P.Con.mu./r.^4.*P.Con.J2.*P.Con.Re.^2;
            n = sqrt(P.Con.mu./a.^3);
            aop = th + w;
            % Forces
            f_r = MJ2/2.*(1-3*sin(inc).^2.*sin(aop).^2);
            f_th = MJ2.*sin(inc).^2.*sin(aop).*cos(aop);
            f_h = MJ2.*sin(inc).*cos(inc).*sin(aop);
            F = reshape([f_r;f_th;f_h],3*P.Con.N_sats,1);
            % GVE Sparsified Matrix
            B11 = vec2sparse(2*a.^2./h.*e.*sin(th),[6,3],[1,1]);
            B12 = vec2sparse(2*a.^2./h.*(1 + e.*cos(th)),[6,3],[1,2]);
            B21 = vec2sparse(p./h.*sin(th),[6,3],[2,1]);
            B22 = vec2sparse(r./h.*(e + 2*cos(th) + e.*cos(th).^2),[6,3],[2,2]);
            B33 = vec2sparse(r./h.*cos(aop),[6,3],[3,3]);
            B43 = vec2sparse(r.*sin(aop)./(h.*sin(inc)),[6,3],[4,3]);
            B51 = vec2sparse(-p./(h.*e).*cos(th),[6,3],[5,1]);
            B52 = vec2sparse(r./(h.*e).*(2+e.*cos(th)).*sin(th),[6,3],[5,2]);
            B53 = vec2sparse(-r./h.*sin(aop).*cos(inc)./sin(inc),[6,3],[5,3]);
            B61 = vec2sparse((p.*cos(th)-2*r.*e)./(n.*a.^2.*e),[6,3],[6,1]);
            B62 = vec2sparse(-(p+r).*sin(th)./(n.*a.^2.*e),[6,3],[6,2]);
            
            B = B11 + B12 + B21 + B22 + B33 + B43 + B51 + B52 + B53 + B61 + B62;
            
            % Kepler Solution
            K = reshape([zeros(5,P.Con.N_sats);n],6*P.Con.N_sats,1);
            
            % Equations of Motion
            dX = B*F + K;
        end
        
%         function dX = dyn_OE_Osc2(P,t,X) 
%             
%         end
        
        
    end
    
    methods(Access = protected, Sealed)
        function propgroups = getPropertyGroups(P)  %#ok<MANU>
            propgroups = matlab.mixin.util.PropertyGroup;
            propgroups(1).Title = 'Propagation Definitions';
            propgroups(1).PropertyList = {'reltol','abstol','Con'};
        end
    end
end