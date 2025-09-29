function dX = ConDynOeOsc(P,t,X)  %#ok<INUSL>
            % rewrite with LPE, Controlled
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
            f  = pi/180*me2ta(Me*180/pi,e); %  conversion M to f
            % Secondary definitions
            p = a.*(1-e.^2);
            h = sqrt(mu*p);
            r = p./(1+e.*cos(f));
            k2 = -3*mu./r.^3.*P.Con.J2.*Re.^2;
            n = sqrt(mu./a.^3);
            eta = sqrt(1-e.^2);
            aol = f + w;
            % Test vs. GVE
            %             MJ2 = -3*P.Con.mu./r.^4.*P.Con.J2.*P.Con.Re.^2;
            %             % Forces
            %             fR = MJ2/2.*(1-3*sin(inc).^2.*sin(aol).^2);
            %             fTh = MJ2.*sin(inc).^2.*sin(aol).*cos(aol);
            %             fH_si = MJ2.*cos(inc).*sin(aol); % sin(inc) taken out to avoid singularity
            %             % Element Rates
            %             da = 2*a.^2./h.*e.*sin(f).*fR +...
            %                 2*a.^2./h.*(1 + e.*cos(f)).*fTh;
            %             de = p./h.*sin(f).*fR +...
            %                 r./h.*(e + 2*cos(f) + e.*cos(f).^2).*fTh;
            %             di = r./h.*cos(aol).*sin(inc).*fH_si;
            %             dO = r.*sin(aol)./(h).*fH_si;
            %             dw = -p./(h.*e).*cos(f).*fR +...
            %                 r./(h.*e).*(2+e.*cos(f)).*sin(f).*fTh +...
            %                 -r./h.*sin(aol).*cos(inc).*fH_si;
            %             dM = (p.*cos(f)-2*r.*e)./(n.*a.^2.*e).*fR +...
            %                 -(p+r).*sin(f)./(n.*a.^2.*e).*fTh;

            dRde = ((e.^2.*cos(f)+cos(f)+2*e).*(3.*sin(inc).^2.*sin(aol).^2-1) ...
                - sin(f).*(2+e.*cos(f)).*(e.*sin(f).*(3*sin(inc).^2.*sin(aol).^2-1) ...
                - 2*sin(inc).^2.*sin(aol).*cos(aol).*(1+e.*cos(f))))./...
                (2*(1+e.*cos(f)).*(1-e.^2));
            % Element Rates
            da = k2.*h.*(e.*sin(f).*(1-3*sin(inc).^2.*sin(aol).^2) +...
                2*sin(inc).^2.*sin(aol).*cos(aol).*(1+e.*cos(f)))./(n.^2.*a.*p.*r);
            de = (1-e.^2).*da./(2.*a.*e) ...
                - k2.*eta.*sin(inc).^2.*sin(aol).*cos(aol)./(n.*a.^2.*e);
            di = k2.*sin(inc).*cos(inc).*cos(aol).*sin(aol)./(n.*a.^2.*eta);
            dO = k2.*cos(inc).*sin(aol).^2./(n.*a.^2.*eta);
            dw = k2.*eta.*dRde./(n.*a.^2.*e) - cos(inc).*dO;
            dM = k2.*(3.*sin(inc).^2.*sin(aol).^2-1)./(n.*a.^2) ...
                - k2.*eta.^2.*dRde./(n.*a.^2.*e);

            dOe = reshape([da;de;di;dO;dw;dM + n],6*P.Con.nSats,1) + ...
                P.Control.ControlOE(t,OE);
            % Equations of Motion
            dX = dOe;
        end