    function [freq0,lpeSpec] = DynOeFourier(P,t,icM,kMax)
            %% Handle Input
            J2 = P.Con.primary.J2;
            Re = P.Con.primary.Re;
            mu = P.Con.primary.mu;

            % handle elements vector
            nT = length(t);
            
            % Propagator now gives Mean elements
            % icM = osc2meNum(icOsc); % change to numerical mean
            % icM(3:end) = icM(3:end)*pi/180;
            % icOsc(3:end) = icOsc(3:end)*pi/180;
            % handle elements vector
            a = icM(1);
            e = icM(2);
            i = icM(3);

            nMo = sqrt(mu/a^3);
            p = a*(1-e^2);
            aop = icM(5);

            b = (1-sqrt(1-e^2))/e;


            %% constant vectors
            m2 = (0:4).';
            m3 = (0:6).';
            m4 = (0:8).';
            m5 = (0:10).';

            a2 = [1,-4*e,4*e^2+2,-4*e,1].';
            a3 = [1,-6*e,12*e^2+3,-8*e^3-12*e,12*e^2+3,-6*e,1].';
            a4 = [1,-8*e,24*e^2+4,-(32*e^3+24*e),(16*e^4+48*e^2+6),...
                -(32*e^3+24*e),24*e^2+4,-8*e,1].';
            a5 = [1,-10*e,40*e^2+5,-(80*e^3+40*e),(80*e^4+120*e^2+10),...
                -(32*e^5+160*e^3+60*e),(80*e^4+120*e^2+10),-(80*e^3+40*e),40*e^2+5,...
                -10*e,1].';

            b1 = [-1,2*e,0,-2*e,1].';
            b2 = [-1,4*e,-4*e^2-1,0,4*e^2+1,-4*e,1].';
            b3 = [-1,6*e,-12*e^2-2,8*e^3+6*e,0,-8*e^3-6*e,12*e^2+2,-6*e,1].';
            b4 = [-1,8*e,-24*e^2-3,32*e^3+16*e,-16*e^4-24*e^2-2,0,16*e^4+24*e^2+2,...
                -32*e^3-16*e,24*e^2+3,-8*e,1].';

            da2de = [0,-4,8*e,-4,0].';
            da3de = [0,-6,24*e,(-24*e^2-12),24*e,-6,0].';
            da4de = [0,-8,48*e,(-96*e^2-24),(64*e^3+96*e),(-96*e^2-24),48*e,-8,0].';
            da5de = [0,-10,80*e,(-240*e^2-40),(320*e^3+240*e),(-160*e^4-480*e^2-60),...
                (320*e^3+240*e),(-240*e^2-40),80*e,-10,0].';

            db1de = [0,2,0,-2,0].';
            db2de = [0,4,-8*e,0,8*e,-4,0].';
            db3de = [0,6,-24*e,24*e^2+6,0,-24*e^2-6,24*e,-6,0].';
            db4de = [0,8,-48*e,96*e^2+16,-64*e^3-48*e,0,64*e^3+48*e,-96*e^2-16,48*e,-8,0].';


            C = [6*(3*sin(i)^2*cos(aop)^2-1)/(1-e^2)^2;
                (9*e^2*sin(i)^2*cos(aop)^2+3*sin(i)^2*(sin(aop)^2-cos(aop)^2)-3*e^2)/2/(1-e^2)^3.5;
                (3*e^3*sin(i)^2*cos(aop)^2+9*e*sin(i)^2*(sin(aop)^2-cos(aop)^2)-e^3)/4/(1-e^2)^4;
                (9*e^2*sin(i)^2*(sin(aop)^2-cos(aop)^2))/8/(1-e^2)^4.5;
                (3*e^3*sin(i)^2*(sin(aop)^2-cos(aop)^2))/16/(1-e^2)^5];

            dCdi_si = [36*cos(aop)^2*cos(i)/(1-e^2)^2; % all pre-divided by sin(i)
                3*cos(i)*(3*e^2*cos(aop)^2-2*cos(aop)^2+1)/(1-e^2)^3.5;
                3*e*cos(i)*(3+(e^2-6)*cos(aop)^2)/2/(1-e^2)^4;
                -9*e^2*cos(i)*(2*cos(aop)^2-1)/4/(1-e^2)^4.5;
                -3*e^3*cos(i)*(2*cos(aop)^2-1)/8/(1-e^2)^5];

            dCdo = [-36*sin(i)^2*sin(aop)*cos(aop)/(1-e^2)^2;
                3*(-3*e^2+2)*sin(i)^2*sin(aop)*cos(aop)/(1-e^2)^3.5;
                -3*e*(e^2-6)*sin(i)^2*sin(aop)*cos(aop)/2/(1-e^2)^4;
                9*e^2*sin(i)^2*sin(aop)*cos(aop)/2/(1-e^2)^4.5;
                3*e^3*sin(i)^2*sin(aop)*cos(aop)/4/(1-e^2)^5];

            dCdo_si = [-36*sin(i)*sin(aop)*cos(aop)/(1-e^2)^2;
                3*(-3*e^2+2)*sin(i)*sin(aop)*cos(aop)/(1-e^2)^3.5;
                -3*e*(e^2-6)*sin(i)*sin(aop)*cos(aop)/2/(1-e^2)^4;
                9*e^2*sin(i)*sin(aop)*cos(aop)/2/(1-e^2)^4.5;
                3*e^3*sin(i)*sin(aop)*cos(aop)/4/(1-e^2)^5];


            dCde = [24*e*(3*sin(i)^2*cos(aop)^2-1)/(1-e^2)^3;
                3*e*(15*e^2*sin(i)^2*cos(aop)^2 +7*sin(i)^2*sin(aop)^2 -sin(i)^2*cos(aop)^2 -5*e^2 -2)/2/(1-e^2)^4.5;
                ((15*e^4-54*e^2-9)*sin(i)^2*cos(aop)^2 + (63*e^2+9)*sin(i)^2*sin(aop)^2 -5*e^4 -3*e^2)/4/(1-e^2)^5;
                -9*e*sin(i)^2*(cos(aop)^2-sin(aop)^2)*(7*e^2+2)/8/(1-e^2)^5.5;
                -3*e^2*sin(i)^2*(cos(aop)^2-sin(aop)^2)*(7*e^2+3)/16/(1-e^2)^6];


            S = -[6*sin(i)^2*sin(aop)*cos(aop)/2/(1-e^2)^3;
                18*e*sin(i)^2*sin(aop)*cos(aop)/4/(1-e^2)^3.5;
                18*e^2*sin(i)^2*sin(aop)*cos(aop)/8/(1-e^2)^4;
                6*e^3*sin(i)^2*sin(aop)*cos(aop)/16/(1-e^2)^4.5];

            dSdi_si = -[6*cos(i)*sin(aop)*cos(aop)/(1-e^2)^3;
                9*e*cos(i)*sin(aop)*cos(aop)/(1-e^2)^3.5;
                9*e^2*cos(i)*sin(aop)*cos(aop)/2/(1-e^2)^4;
                3*e^3*cos(i)*sin(aop)*cos(aop)/4/(1-e^2)^4.5];

            dSdo = -[6*sin(i)^2*(2*cos(aop)^2-1)/2/(1-e^2)^3;
                18*e*sin(i)^2*(2*cos(aop)^2-1)/4/(1-e^2)^3.5;
                18*e^2*sin(i)^2*(2*cos(aop)^2-1)/8/(1-e^2)^4;
                6*e^3*sin(i)^2*(2*cos(aop)^2-1)/16/(1-e^2)^4.5];

            dSdo_si = -[6*sin(i)*(2*cos(aop)^2-1)/2/(1-e^2)^3;
                18*e*sin(i)*(2*cos(aop)^2-1)/4/(1-e^2)^3.5;
                18*e^2*sin(i)*(2*cos(aop)^2-1)/8/(1-e^2)^4;
                6*e^3*sin(i)*(2*cos(aop)^2-1)/16/(1-e^2)^4.5];

            dSde = [-18*e*sin(i)^2*sin(aop)*cos(aop)/(1-e^2)^4;
                -9*sin(i)^2*sin(aop)*cos(aop)*(6*e^2+1)/2/(1-e^2)^4.5;
                -9*e*sin(i)^2*sin(aop)*cos(aop)*(3*e^2+1)/2/(1-e^2)^5;
                -9*e^2*sin(i)^2*sin(aop)*cos(aop)*(2*e^2+1)/8/(1-e^2)^5.5];

            AkM = nan(5,kMax);
            Ak_eM = nan(5,kMax);
            Akde_eM = nan(5,kMax);

            BkM = nan(4,kMax);
            Bk_eM = nan(4,kMax);
            Bkde_eM = nan(4,kMax);

            k = 1;

            while k <= kMax
                n = 0;

                g2 = b.^abs(m2+n+k-2);
                g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
                g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
                    3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
                    3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
                g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
                    e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
                    5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
                    5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);

                dg2deXe = abs(m2+n+k-2).*b.^(abs(m2+n+k-2))/sqrt(1-e^2);
                dg3deXe = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^(abs(m3+n+k-3)) + ...
                    e/sqrt(1-e^2)*b.^(abs(m3+n+k-2)))/sqrt(1-e^2) + ...
                    e*b.^abs(m3+n+k-2)/(1-e^2)^(3/2);
                dg4deXe = abs(m4+n+k-2).*(3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2) +...
                    abs(m4+n+k-3).*(3*e/2/sqrt(1-e^2)*b.^abs(m4+n+k-3) +...
                    abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))/sqrt(1-e^2) +...
                    3/2*e*abs(m4+n+k-2)/(1-e^2)^(3/2).*b.^abs(m4+n+k-3) +...
                    3*e^2/(1-e^2)^2*b.^abs(m4+n+k-2);
                dg5deXe = abs(m5+n+k-2).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(m5+n+k-2) + ...
                    abs(m5+n+k-3).*(5*e^2/2/(1-e^2)*b.^abs(m5+n+k-3) + ...
                    abs(m5+n+k-4).*(e/sqrt(1-e^2)*b.^abs(m5+n+k-4) + ...
                    abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))/sqrt(1-e^2) + ...
                    abs(m5+n+k-3).*abs(m5+n+k-2)*e/(1-e^2)^(3/2).*b.^abs(m5+n+k-4) +...
                    5*abs(m5+n+k-2)*e^2/(1-e^2)^2.*b.^abs(m5+n+k-3) +...
                    15*e^3/2/(1-e^2)^(5/2)*b.^abs(m5+n+k-2);

                Jk = besselj(k,k*e);
                % Jk/e nonsingular
                Jk_e = 0.5*(besselj(k+1,k*e) + besselj(k-1,k*e));
                % dJkde/e nonsingular
                if k~=1
                    % use expression with no /e
                    dJkde_e = k^2/4/(k^2-1)*(2*Jk + (k+1)*besselj(k-2,k*e) - (k-1)*besselj(k+2,k*e));
                else
                    % elimination of /e not possible
                    dJkde_e = 0.5*k/e*(besselj(k-1,k*e) - besselj(k+1,k*e));
                end

                Jn = besselj(n,-k*e);
                % Jn/e nonsingular
                if n~=0
                    Jn_e = -k/2/n*(besselj(n+1,-k*e)+besselj(n-1,-k*e));
                else
                    Jn_e = Jn/e;
                end
                % dJnde/e nonsingular
                if abs(n)~=1
                    dJnde_e = k^2/4/(n^2-1)*(2*Jn + (n+1)*besselj(n-2,-k*e) - (n-1)*besselj(n+2,-k*e));
                else
                    dJnde_e = -0.5*k/e*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
                end
                % Jn/e^2 nonsingular
                if abs(n)>1
                    Jn_e2 = k^2/4/n/(n^2-1)*(2*n*Jn + (n+1)*besselj(n-2,-k*e) + (n-1)*besselj(n+2,-k*e));
                else
                    Jn_e2 = Jn/e^2;
                end

                Ak = [Jk; Jn*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                Ak_e = [Jk_e; Jn_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                Akde_e = [dJkde_e; dJnde_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] +...
                    Jn_e*[da2de.'*g2; da3de.'*g3; da4de.'*g4; da5de.'*g5] +...
                    Jn_e2*[a2.'*dg2deXe; a3.'*dg3deXe; a4.'*dg4deXe; a5.'*dg5deXe]];

                Bk = Jn*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                Bk_e = Jn_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                Bkde_e = dJnde_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] +...
                    Jn_e*[db1de.'*g2; db2de.'*g3; db3de.'*g4; db4de.'*g5] +...
                    Jn_e2*[b1.'*dg2deXe; b2.'*dg3deXe; b3.'*dg4deXe; b4.'*dg5deXe];

                n = 1;
                while n <= k + 3
                    % positive n
                    g2 = b.^abs(m2+n+k-2);
                    g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
                    g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
                        3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
                        3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
                    g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
                        e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
                        5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
                        5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);

                    dg2deXe = abs(m2+n+k-2).*b.^(abs(m2+n+k-2))/sqrt(1-e^2);
                    dg3deXe = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^(abs(m3+n+k-3)) + ...
                        e/sqrt(1-e^2)*b.^(abs(m3+n+k-2)))/sqrt(1-e^2) + ...
                        e*b.^abs(m3+n+k-2)/(1-e^2)^(3/2);
                    dg4deXe = abs(m4+n+k-2).*(3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2) +...
                        abs(m4+n+k-3).*(3*e/2/sqrt(1-e^2)*b.^abs(m4+n+k-3) +...
                        abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))/sqrt(1-e^2) +...
                        3/2*e*abs(m4+n+k-2)/(1-e^2)^(3/2).*b.^abs(m4+n+k-3) +...
                        3*e^2/(1-e^2)^2*b.^abs(m4+n+k-2);
                    dg5deXe = abs(m5+n+k-2).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(m5+n+k-2) + ...
                        abs(m5+n+k-3).*(5*e^2/2/(1-e^2)*b.^abs(m5+n+k-3) + ...
                        abs(m5+n+k-4).*(e/sqrt(1-e^2)*b.^abs(m5+n+k-4) + ...
                        abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))/sqrt(1-e^2) + ...
                        abs(m5+n+k-3).*abs(m5+n+k-2)*e/(1-e^2)^(3/2).*b.^abs(m5+n+k-4) +...
                        5*abs(m5+n+k-2)*e^2/(1-e^2)^2.*b.^abs(m5+n+k-3) +...
                        15*e^3/2/(1-e^2)^(5/2)*b.^abs(m5+n+k-2);

                    Jn = besselj(n,-k*e);
                    % Jn/e nonsingular
                    if n~=0
                        Jn_e = -k/2/n*(besselj(n+1,-k*e)+besselj(n-1,-k*e));
                    else
                        Jn_e = Jn/e;
                    end
                    % dJnde/e nonsingular
                    if abs(n)~=1
                        dJnde_e = k^2/4/(n^2-1)*(2*Jn + (n+1)*besselj(n-2,-k*e) - (n-1)*besselj(n+2,-k*e));
                    else
                        dJnde_e = -0.5*k/e*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
                    end
                    % Jn/e^2 nonsingular
                    if abs(n)>1
                        Jn_e2 = k^2/4/n/(n^2-1)*(2*n*Jn + (n+1)*besselj(n-2,-k*e) + (n-1)*besselj(n+2,-k*e));
                    else
                        Jn_e2 = Jn/e^2;
                    end

                    dAk = [0; Jn*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                    dAk_e = [0; Jn_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                    dAkde_e = [0; dJnde_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] +...
                        Jn_e*[da2de.'*g2; da3de.'*g3; da4de.'*g4; da5de.'*g5] +...
                        Jn_e2*[a2.'*dg2deXe; a3.'*dg3deXe; a4.'*dg4deXe; a5.'*dg5deXe]];

                    dBk = Jn*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBk_e = Jn_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBkde_e = dJnde_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] +...
                        Jn_e*[db1de.'*g2; db2de.'*g3; db3de.'*g4; db4de.'*g5] +...
                        Jn_e2*[b1.'*dg2deXe; b2.'*dg3deXe; b3.'*dg4deXe; b4.'*dg5deXe];


                    % negative n
                    n = -n;
                    g2 = b.^abs(m2+n+k-2);
                    g3 = abs(m3+n+k-2).*b.^abs(m3+n+k-3) + e/sqrt(1-e^2)*b.^abs(m3+n+k-2);
                    g4 = abs(m4+n+k-3).*abs(m4+n+k-2)/2.*b.^abs(m4+n+k-4) + ...
                        3*e*abs(m4+n+k-2)/2/sqrt(1-e^2).*b.^abs(m4+n+k-3) + ...
                        3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2);
                    g5 = abs(m5+n+k-4).*abs(m5+n+k-3).*abs(m5+n+k-2)/6.*b.^abs(m5+n+k-5) + ...
                        e*abs(m5+n+k-3).*abs(m5+n+k-2)/sqrt(1-e^2).*b.^abs(m5+n+k-4) + ...
                        5*e^2*abs(m5+n+k-2)/2/(1-e^2).*b.^abs(m5+n+k-3) + ...
                        5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5+n+k-2);

                    dg2deXe = abs(m2+n+k-2).*b.^(abs(m2+n+k-2))/sqrt(1-e^2);
                    dg3deXe = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^(abs(m3+n+k-3)) + ...
                        e/sqrt(1-e^2)*b.^(abs(m3+n+k-2)))/sqrt(1-e^2) + ...
                        e*b.^abs(m3+n+k-2)/(1-e^2)^(3/2);
                    dg4deXe = abs(m4+n+k-2).*(3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2) +...
                        abs(m4+n+k-3).*(3*e/2/sqrt(1-e^2)*b.^abs(m4+n+k-3) +...
                        abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))/sqrt(1-e^2) +...
                        3/2*e*abs(m4+n+k-2)/(1-e^2)^(3/2).*b.^abs(m4+n+k-3) +...
                        3*e^2/(1-e^2)^2*b.^abs(m4+n+k-2);
                    dg5deXe = abs(m5+n+k-2).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(m5+n+k-2) + ...
                        abs(m5+n+k-3).*(5*e^2/2/(1-e^2)*b.^abs(m5+n+k-3) + ...
                        abs(m5+n+k-4).*(e/sqrt(1-e^2)*b.^abs(m5+n+k-4) + ...
                        abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))/sqrt(1-e^2) + ...
                        abs(m5+n+k-3).*abs(m5+n+k-2)*e/(1-e^2)^(3/2).*b.^abs(m5+n+k-4) +...
                        5*abs(m5+n+k-2)*e^2/(1-e^2)^2.*b.^abs(m5+n+k-3) +...
                        15*e^3/2/(1-e^2)^(5/2)*b.^abs(m5+n+k-2);

                    Jn = besselj(n,-k*e);
                    % Jn/e nonsingular
                    if n~=0
                        Jn_e = -k/2/n*(besselj(n+1,-k*e)+besselj(n-1,-k*e));
                    else
                        Jn_e = Jn/e;
                    end
                    % dJnde/e nonsingular
                    if abs(n)~=1
                        dJnde_e = k^2/4/(n^2-1)*(2*Jn + (n+1)*besselj(n-2,-k*e) - (n-1)*besselj(n+2,-k*e));
                    else
                        dJnde_e = -0.5*k/e*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
                    end
                    % Jn/e^2 nonsingular
                    if abs(n)>1
                        Jn_e2 = k^2/4/n/(n^2-1)*(2*n*Jn + (n+1)*besselj(n-2,-k*e) + (n-1)*besselj(n+2,-k*e));
                    else
                        Jn_e2 = Jn/e^2;
                    end

                    dAk = dAk + [0; Jn*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                    dAk_e = dAk_e + [0; Jn_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]];
                    dAkde_e = dAkde_e + [0; dJnde_e*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] +...
                        Jn_e*[da2de.'*g2; da3de.'*g3; da4de.'*g4; da5de.'*g5] +...
                        Jn_e2*[a2.'*dg2deXe; a3.'*dg3deXe; a4.'*dg4deXe; a5.'*dg5deXe]];

                    dBk = dBk + Jn*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBk_e = dBk_e + Jn_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
                    dBkde_e = dBkde_e + dJnde_e*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] +...
                        Jn_e*[db1de.'*g2; db2de.'*g3; db3de.'*g4; db4de.'*g5] +...
                        Jn_e2*[b1.'*dg2deXe; b2.'*dg3deXe; b3.'*dg4deXe; b4.'*dg5deXe];

                    Ak = Ak + dAk;
                    Ak_e = Ak_e + dAk_e;
                    Akde_e = Akde_e + dAkde_e;

                    Bk = Bk + dBk;
                    Bk_e = Bk_e + dBk_e;
                    Bkde_e = Bkde_e + dBkde_e;

                    n = abs(n); % return n to positive value
                    n = n+1;
                end

                AkM(:,k) = Ak;
                Ak_eM(:,k) = Ak_e;
                Akde_eM(:,k) = Akde_e;

                BkM(:,k) = Bk;
                Bk_eM(:,k) = Bk_e;
                Bkde_eM(:,k) = Bkde_e;

                k = k+1;
            end


            eta = sqrt(1-e^2);
            % constant potential values
            R = -nMo^2*Re^2/2; % common factor
            Rsma = -nMo*J2*Re^2/a;
            Recc = -nMo*eta*J2*Re^2/2/a^2;
            Rinc = -nMo*J2*Re^2*cos(i)/2/a^2/eta;
            Rran = -nMo*J2*Re^2/2/a^2/eta;
            Raop = -nMo*J2*Re^2/2/a^2;
            Rman = Raop;

            % Freq 0 elements without common factor
            ran0 = 3*cos(i)/eta^3;
            aop0 = -1.5*(5*cos(i)^2-1)/eta^4;
            man0 = -1.5*(3*cos(i)^2-1)/eta^3;

            % Calculate Spectrum of Elements
            k = 1:kMax;
            lpeSpec = nan(12,kMax);

            lpeSpec(1:2,:) = Rsma*[S.'*(BkM.*k);
                -C.'*(AkM.*k)];
            lpeSpec(3:4,:) = Recc*[eta*S.'*(Bk_eM.*k) - dCdo.'*Ak_eM;
                -eta*C.'*(Ak_eM.*k) - dSdo.'*Bk_eM];
            lpeSpec(5:6,:) = Rinc*[dCdo_si.'*AkM;
                dSdo_si.'*BkM];
            lpeSpec(7:8,:) = Rran*[dCdi_si.'*AkM;
                dSdi_si.'*BkM];
            lpeSpec(9:10,:) = Raop*[eta*(dCde.'*Ak_eM + C.'*Akde_eM) - cos(i)/eta*dCdi_si.'*AkM;
                eta*(dSde.'*Bk_eM + S.'*Bkde_eM) - cos(i)/eta*dSdi_si.'*BkM];
            lpeSpec(11:12,:) = Rman*[-eta^2*(dCde.'*Ak_eM + C.'*Akde_eM) + 6*C.'*AkM;
                -eta^2*(dSde.'*Bk_eM + S.'*Bkde_eM) + 6*S.'*BkM];
            %% Zero Frequency Components
            freq0 = [0; 0; 0; Rran*ran0; Raop*aop0; Rman*man0];
        end