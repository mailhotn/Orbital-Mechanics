function dX = DynOePns(P,t,X)  %#ok<INUSL>
            % Dynamics for Nonsingular polar-nodal propagation
            mu = 1;
            J2 = P.Con.J2;
            Re = 1;
            PNS  = reshape(X,6,P.Con.nSats);
            % Element vectors (angles already in radians)
            r   = PNS(1,:); % radius
            th  = PNS(2,:); % mod. argument of latitude
            n   = PNS(3,:); % mod. right ascension
            R   = PNS(4,:); % Radial velocity
            TH  = PNS(5,:); % specific angular momentum
            N   = PNS(6,:); % mod. angular momentum in z

            % remove singularities n/N = nan iff N=0 i.e. orbit is equatorial
            nbyN = n./N;
            nbyN(isnan(nbyN)) = 0;
            k2 = J2*Re^2/2;

            % Keplerian Motion
            dKep = [R;
                TH./r.^2;
                zeros(1,P.Con.nSats);
                TH.^2./r.^3-mu./r.^2;
                zeros(1,P.Con.nSats);
                zeros(1,P.Con.nSats)];
            den = (2*r.^3.*TH.^2);
            dJ2 = 3*mu*k2*[zeros(1,P.Con.nSats);
                (N.^4-2*TH.*N.^2).*sin(th+nbyN).^2./(den.*TH);
                -sin(th+nbyN).*((2*N.^3-4*TH.*N).*sin(th+nbyN) -cos(th+nbyN).*n.*(N.^2-4*TH))./(den);
                -3*((N.^4-4*TH.*N.^2).*sin(th+nbyN).^2 +4*TH.^2/3)./(2*r.*den);
                (N.^4-4*TH.*N.^2).*sin(th+nbyN).*cos(th+nbyN)./(den);
                (N.^3-4*TH.*N).*sin(th+nbyN).*cos(th+nbyN)./(den)];
            dOe = reshape(dKep + dJ2,6*P.Con.nSats,1);
            % Equations of Motion
            dX = dOe;
        end