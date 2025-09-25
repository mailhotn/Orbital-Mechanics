 function dX = DynOeOsc(P,t,X)  %#ok<*INUSD>
            OE  = reshape(X,6,P.Con.nSats);
            % Normalized Parameters
            mu = 1;
            Re = 1;
            % Element vectors (angles already in radians)
            a   = OE(1,:);
            e   = OE(2,:);
            inc = OE(3,:);
            w   = OE(5,:);
            Me  = OE(6,:);
            th  = pi/180*me2ta(Me*180/pi,e);
            % Secondary definitions
            p = a.*(1-e.^2);
            h = sqrt(mu*p);
            r = p./(1+e.*cos(th));
            MJ2 = -3*mu./r.^4.*P.Con.J2.*Re.^2;
            n = sqrt(mu./a.^3);
            aol = th + w;
            % Forces
            fR = MJ2/2.*(1-3*sin(inc).^2.*sin(aol).^2);
            fTh = MJ2.*sin(inc).^2.*sin(aol).*cos(aol);
            fH_si = MJ2.*cos(inc).*sin(aol); % sin(inc) taken out to avoid singularity
            % Element Rates
            da = 2*a.^2./h.*e.*sin(th).*fR +...
                2*a.^2./h.*(1 + e.*cos(th)).*fTh;
            de = p./h.*sin(th).*fR +...
                r./h.*(e + 2*cos(th) + e.*cos(th).^2).*fTh;
            di = r./h.*cos(aol).*sin(inc).*fH_si;
            dO = r.*sin(aol)./(h).*fH_si;
            dw = -p./(h.*e).*cos(th).*fR +...
                r./(h.*e).*(2+e.*cos(th)).*sin(th).*fTh +...
                -r./h.*sin(aol).*cos(inc).*fH_si;
            dM = (p.*cos(th)-2*r.*e)./(n.*a.^2.*e).*fR +...
                -(p+r).*sin(th)./(n.*a.^2.*e).*fTh;

            dOe = reshape([da;de;di;dO;dw;dM + n],6*P.Con.nSats,1);
            % Equations of Motion
            dX = dOe;
        end