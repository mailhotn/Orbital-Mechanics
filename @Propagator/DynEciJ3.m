function dX = DynEciJ3(P,t,X) %#ok<INUSL>
            % move stuff around
            order = 6*P.Con.nSats;
            X2 = reshape(X,[6,P.Con.nSats]);
            Rv = [eye(3),zeros(3);
                zeros(3,6)]*X2;
            % get vector of R magnitudes
            r = sqrt(dot(Rv,Rv,1));
            r = reshape(repmat(r,6,1),order,1);
            % move more stuff around
            R2 = repmat([1 1 1 0 0 0].',P.Con.nSats,1).*X;
            V2 = repmat([0 0 0 1 1 1].',P.Con.nSats,1).*X;
            Z  = repmat([0 0 1 0 0 0].',P.Con.nSats,1).*R2;
            Z2 = [zeros(3,2),ones(3,1),zeros(3);zeros(3,6)]*X2;
            Z2 = reshape(Z2,[order,1]);
            % equation of motion
            f_J2 = -P.Con.mu*P.Con.J2*P.Con.Re^2./r.^4.*...
                (3*Z./r + (-7.5*(Z2./r).^2 + 1.5).*R2./r);

            f_J3 = -5/2*P.Con.mu*P.Con.J3*P.Con.Re^3./r.^7.*...
                ((3*Z2 - 7*Z2.^2./r.^2).*R2 + ...
                (3 - 3/5*r.^2./(Z2.^2)).*Z.^2);

            dX = circshift(V2,-3) + ...
                circshift(-P.Con.mu*R2./r.^3 + f_J2 + f_J3,3);
        end