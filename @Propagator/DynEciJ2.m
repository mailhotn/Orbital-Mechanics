function dX = DynEciJ2(P,t,X) %#ok<INUSL>
            % Normalized Parameters
            mu = 1;
            Re = 1;
            % move stuff around
            order = 6*P.Con.nSats;
            X2 = reshape(X,[6,P.Con.nSats]);
            Rv = [eye(3),zeros(3);
                zeros(3,6)]*X2;
            % get vector of R magnitudes
            %             r = sqrt(dot(Rv,Rv,1));
            r = vecnorm(Rv);
            r = reshape(repmat(r,6,1),order,1);
            % move more stuff around
            R2 = reshape(Rv,order,1);
            V2 = reshape([zeros(3,6);zeros(3),eye(3)]*X2,order,1);
            Z  = repmat([0 0 1 0 0 0].',P.Con.nSats,1).*R2;
            Z2 = reshape([zeros(3,2),ones(3,1),zeros(3);zeros(3,6)]*X2,order,1);
            % equation of motion
            f_J2 = -mu*P.Con.J2*Re^2./r.^4.*...
                (3*Z./r + (-7.5*(Z2./r).^2 + 1.5).*R2./r);
            dX = circshift(V2,-3) + circshift(-mu*R2./r.^3 + f_J2,3);
        end