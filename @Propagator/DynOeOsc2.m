function dX = DynOeOsc2(P,t,X)  %#ok<INUSL>
            OE  = reshape(X,6,P.Con.nSats);
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
            aol = th + w;
            % Forces
            f_r = MJ2/2.*(1-3*sin(inc).^2.*sin(aol).^2);
            f_th = MJ2.*sin(inc).^2.*sin(aol).*cos(aol);
            f_h = MJ2.*sin(inc).*cos(inc).*sin(aol);
            F = reshape([f_r;f_th;f_h],3*P.Con.nSats,1);
            % GVE Sparsified Matrix
            B11 = vec2sparse(2*a.^2./h.*e.*sin(th),[6,3],[1,1]);
            B12 = vec2sparse(2*a.^2./h.*(1 + e.*cos(th)),[6,3],[1,2]);
            B21 = vec2sparse(p./h.*sin(th),[6,3],[2,1]);
            B22 = vec2sparse(r./h.*(e + 2*cos(th) + e.*cos(th).^2),[6,3],[2,2]);
            B33 = vec2sparse(r./h.*cos(aol),[6,3],[3,3]);
            B43 = vec2sparse(r.*sin(aol)./(h.*sin(inc)),[6,3],[4,3]);
            B51 = vec2sparse(-p./(h.*e).*cos(th),[6,3],[5,1]);
            B52 = vec2sparse(r./(h.*e).*(2+e.*cos(th)).*sin(th),[6,3],[5,2]);
            B53 = vec2sparse(-r./h.*sin(aol).*cos(inc)./sin(inc),[6,3],[5,3]);
            B61 = vec2sparse((p.*cos(th)-2*r.*e)./(n.*a.^2.*e),[6,3],[6,1]);
            B62 = vec2sparse(-(p+r).*sin(th)./(n.*a.^2.*e),[6,3],[6,2]);

            B = B11 + B12 + B21 + B22 + B33 + B43 + B51 + B52 + B53 + B61 + B62;

            % Kepler Solution
            K = reshape([zeros(5,P.Con.nSats);n],6*P.Con.nSats,1);

            % Equations of Motion
            dX = B*F + K;
        end