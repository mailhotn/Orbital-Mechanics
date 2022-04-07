function [dOe,dJ2] = LpeJ2Fourier(t,oe,kMax,primary)

if nargin < 4
    primary = earth();
end
% handle elements vector
a = oe(1);
e = oe(2);
i = oe(3);
raan = oe(4);
aop = oe(5);
M = oe(6);
b = (1-sqrt(1-e^2))/e;

% constant potential values
R = -primary.mu*primary.J2*primary.Re^2/2/a^3; % common factor
dRda = -3*R/a;

R0 = -(3*cos(i)^2-1)/2/(1-e^2)^(3/2); % 0 freq element
dR0di = 3*cos(i)*sin(i)/(1-e^2)^(3/2);
dR0de = -3*e*(3*cos(i)^2-1)/2/(1-e^2)^(5/2); %
dR0do = 0;
dR0dl = 0;

% constant vectors
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


C1xA1 = 6*(3*sin(i)^2*cos(aop)^2-1)/(1-e^2)^2;

dC1xA1di = 36*sin(i)*cos(i)*cos(aop)^2/(1-e^2)^2;
dC1xA1do = -36*sin(i)^2*sin(aop)*cos(aop)/(1-e^2)^2;
dC1xA1de =  24*e*(3*sin(i)^2*cos(aop)^2-1)/(1-e^2)^3;

C = [(9*e^2*sin(i)^2*cos(aop)^2+3*sin(i)^2*(sin(aop)^2-cos(aop)^2)-3*e^2)/2/(1-e^2)^(7/2);
    (3*e^3*sin(i)^2*cos(aop)^2+9*e*sin(i)^2*(sin(aop)^2-cos(aop)^2)-e^3)/4/(1-e^2)^4;
    (9*e^2*sin(i)^2*(sin(aop)^2-cos(aop)^2))/8/(1-e^2)^(9/2);
    (3*e^3*sin(i)^2*(sin(aop)^2-cos(aop)^2))/16/(1-e^2)^5];

dCdi = [3*(3*e^2*sin(i)*cos(i)*cos(aop)^2+sin(i)*cos(i)*(sin(aop)^2-cos(aop)^2))/(1-e^2)^(7/2);
    3*(e^3*sin(i)*cos(i)*cos(aop)^2+3*e*sin(i)*cos(i)*(sin(aop)^2-cos(aop)^2))/2/(1-e^2)^4;
    9*e^2*sin(i)*cos(i)*(sin(aop)^2-cos(aop)^2)/4/(1-e^2)^(9/2);
    3*e^3*sin(i)*cos(i)*(sin(aop)^2-cos(aop)^2)/8/(1-e^2)^5];

dCdo = [3*(-3*e^2+2)*sin(i)^2*sin(aop)*cos(aop)/(1-e^2)^(7/2);
    -3*e*(e^2-6)*sin(i)^2*sin(aop)*cos(aop)/2/(1-e^2)^4;
    9*e^2*sin(i)^2*sin(aop)*cos(aop)/2/(1-e^2)^(9/2);
    3*e^3*sin(i)^2*sin(aop)*cos(aop)/4/(1-e^2)^5];

dCde = [3*e*(15*e^2*sin(i)^2*cos(aop)^2 +7*sin(i)^2*sin(aop)^2 -sin(i)^2*cos(aop)^2 -5*e^2 -2)/2/(1-e^2)^(9/2);
    ((15*e^4-54*e^2-9)*sin(i)^2*cos(aop)^2 + (63*e^2+9)*sin(i)^2*sin(aop)^2 -5*e^4 -3*e^2)/4/(1-e^2)^5;
    -9*e*sin(i)^2*(cos(aop)^2-sin(aop)^2)*(7*e^2+2)/8/(1-e^2)^(11/2);
    -3*e^2*sin(i)^2*(cos(aop)^2-sin(aop)^2)*(7*e^2+3)/16/(1-e^2)^6];


S = -[6*sin(i)^2*sin(aop)*cos(aop)/2/(1-e^2)^3;
    18*e*sin(i)^2*sin(aop)*cos(aop)/4/(1-e^2)^(7/2);
    18*e^2*sin(i)^2*sin(aop)*cos(aop)/8/(1-e^2)^4;
    6*e^3*sin(i)^2*sin(aop)*cos(aop)/16/(1-e^2)^(9/2)];

dSdi = -[12*sin(i)*cos(i)*sin(aop)*cos(aop)/2/(1-e^2)^3;
    36*e*sin(i)*cos(i)*sin(aop)*cos(aop)/4/(1-e^2)^(7/2);
    36*e^2*sin(i)*cos(i)*sin(aop)*cos(aop)/8/(1-e^2)^4;
    12*e^3*sin(i)*cos(i)*sin(aop)*cos(aop)/16/(1-e^2)^(9/2)];

dSdo = -[6*sin(i)^2*(2*cos(aop)^2-1)/2/(1-e^2)^3;
    18*e*sin(i)^2*(2*cos(aop)^2-1)/4/(1-e^2)^(7/2);
    18*e^2*sin(i)^2*(2*cos(aop)^2-1)/8/(1-e^2)^4;
    6*e^3*sin(i)^2*(2*cos(aop)^2-1)/16/(1-e^2)^(9/2)];

dSde = [-18*e*sin(i)^2*sin(aop)*cos(aop)/(1-e^2)^4;
    -9*sin(i)^2*sin(aop)*cos(aop)*(6*e^2+1)/2/(1-e^2)^(9/2);
    -9*e*sin(i)^2*sin(aop)*cos(aop)*(3*e^2+1)/2/(1-e^2)^5;
    -9*e^2*sin(i)^2*sin(aop)*cos(aop)*(2*e^2+1)/8/(1-e^2)^(11/2)];

dJ2daFreq = dRda*[[R0;0],zeros(2,kMax)];
dJ2deFreq = R*[[dR0de;0],zeros(2,kMax)];
dJ2diFreq = R*[[dR0di;0],zeros(2,kMax)];
dJ2doFreq = R*[[dR0do;0],zeros(2,kMax)];
dJ2dlFreq = R*[[dR0dl;0],zeros(2,kMax)];

k = 1;

while k < kMax
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
    
    dg2de = abs(m2+n+k-2).*b.^(abs(m2+n+k-2))/e/sqrt(1-e^2);
    dg3de = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^(abs(m3+n+k-3)) + ...
        e/sqrt(1-e^2)*b.^(abs(m3+n+k-2)))/e/sqrt(1-e^2) + ...
        b.^abs(m3+n+k-2)/(1-e^2)^(3/2);
    dg4de = abs(m4+n+k-2).*(3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2) +...
        abs(m4+n+k-3).*(3*e/2/sqrt(1-e^2)*b.^abs(m4+n+k-3) +...
        abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))/e/sqrt(1-e^2) +...
        3/2*abs(m4+n+k-2)/(1-e^2)^(3/2).*b.^abs(m4+n+k-3) +...
        3*e/(1-e^2)^2*b.^abs(m4+n+k-2);
    dg5de = abs(m5+n+k-2).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(m5+n+k-2) + ...
        abs(m5+n+k-3).*(5*e^2/2/(1-e^2)*b.^abs(m5+n+k-3) + ...
        abs(m5+n+k-4).*(e/sqrt(1-e^2)*b.^abs(m5+n+k-4) + ...
        abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))/e/sqrt(1-e^2) + ...
        abs(m5+n+k-3).*abs(m5+n+k-2)/(1-e^2)^(3/2).*b.^abs(m5+n+k-4) +...
        5*abs(m5+n+k-2)*e/(1-e^2)^2.*b.^abs(m5+n+k-3) +...
        15*e^2/2/(1-e^2)^(5/2)*b.^abs(m5+n+k-2);
    
    Jk = besselj(k,k*e);
    dJkde = 0.5*(besselj(k-1,k*e) - besselj(k+1,k*e));
    Jn = besselj(n,-k*e);
    dJnde = 0.5*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
    
    Akda = C1xA1*Jk +...
        Jn*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
    Akde = (C1xA1*k*dJkde + dC1xA1de*Jk) + ...
        -k*dJnde*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]+...
        Jn*(dCde.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] + ...
        C.'*[da2de.'*g2+a2.'*dg2de; da3de.'*g3+a3.'*dg3de; da4de.'*g4+a4.'*dg4de; da5de.'*g5+a5.'*dg5de]);
    Akdi = dC1xA1di*Jk +...
        Jn*dCdi.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
    Akdo = dC1xA1do*Jk +...
        Jn*dCdo.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
    Akdl = Akda;
    
    Bkda = Jn*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
    Bkde = -k*dJnde*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5]+...
        Jn*(dSde.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] + ...
        S.'*[db1de.'*g2+b1.'*dg2de; db2de.'*g3+b2.'*dg3de; db3de.'*g4+b3.'*dg4de; db4de.'*g5+b4.'*dg5de]);
    Bkdi = Jn*dSdi.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
    Bkdo = Jn*dSdo.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
    Bkdl = Bkda;
    
    n = 1;
    while n <= k + 8
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
        
        dg2de = abs(m2+n+k-2).*b.^(abs(m2+n+k-2))/e/sqrt(1-e^2);
        dg3de = abs(m3+n+k-2).*(abs(m3+n+k-3).*b.^(abs(m3+n+k-3)) + ...
            e/sqrt(1-e^2)*b.^(abs(m3+n+k-2)))/e/sqrt(1-e^2) + ...
            b.^abs(m3+n+k-2)/(1-e^2)^(3/2);
        dg4de = abs(m4+n+k-2).*(3*e^2/2/(1-e^2)*b.^abs(m4+n+k-2) +...
            abs(m4+n+k-3).*(3*e/2/sqrt(1-e^2)*b.^abs(m4+n+k-3) +...
            abs(m4+n+k-4)/2.*b.^abs(m4+n+k-4)))/e/sqrt(1-e^2) +...
            3/2*abs(m4+n+k-2)/(1-e^2)^(3/2).*b.^abs(m4+n+k-3) +...
            3*e/(1-e^2)^2*b.^abs(m4+n+k-2);
        dg5de = abs(m5+n+k-2).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(m5+n+k-2) + ...
            abs(m5+n+k-3).*(5*e^2/2/(1-e^2)*b.^abs(m5+n+k-3) + ...
            abs(m5+n+k-4).*(e/sqrt(1-e^2)*b.^abs(m5+n+k-4) + ...
            abs(m5+n+k-5)/6.*b.^abs(m5+n+k-5))))/e/sqrt(1-e^2) + ...
            abs(m5+n+k-3).*abs(m5+n+k-2)/(1-e^2)^(3/2).*b.^abs(m5+n+k-4) +...
            5*abs(m5+n+k-2)*e/(1-e^2)^2.*b.^abs(m5+n+k-3) +...
            15*e^2/2/(1-e^2)^(5/2)*b.^abs(m5+n+k-2);
        
        Jn = besselj(n,-k*e);
        dJnde = 0.5*(besselj(n-1,-k*e) - besselj(n+1,-k*e));
        
        dAkda = Jn*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dAkde = -k*dJnde*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]+...
            Jn*(dCde.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] + ...
            C.'*[da2de.'*g2+a2.'*dg2de; da3de.'*g3+a3.'*dg3de; da4de.'*g4+a4.'*dg4de; da5de.'*g5+a5.'*dg5de]);
        dAkdi = Jn*dCdi.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dAkdo = Jn*dCdo.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dAkdl = dAkda;
        
        dBkda = Jn*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        dBkde = -k*dJnde*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5]+...
            Jn*(dSde.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] + ...
            S.'*[db1de.'*g2+b1.'*dg2de; db2de.'*g3+b2.'*dg3de; db3de.'*g4+b3.'*dg4de; db4de.'*g5+b4.'*dg5de]);
        dBkdi = Jn*dSdi.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        dBkdo = Jn*dSdo.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        dBkdl = dBkda;
        
        % negative n
        g2 = b.^abs(m2-n+k-2);
        g3 = abs(m3-n+k-2).*b.^abs(m3-n+k-3) + e/sqrt(1-e^2)*b.^abs(m3-n+k-2);
        g4 = abs(m4-n+k-3).*abs(m4-n+k-2)/2.*b.^abs(m4-n+k-4) + ...
            3*e*abs(m4-n+k-2)/2/sqrt(1-e^2).*b.^abs(m4-n+k-3) + ...
            3*e^2/2/(1-e^2)*b.^abs(m4-n+k-2);
        g5 = abs(m5-n+k-4).*abs(m5-n+k-3).*abs(m5-n+k-2)/6.*b.^abs(m5-n+k-5) + ...
            e*abs(m5-n+k-3).*abs(m5-n+k-2)/sqrt(1-e^2).*b.^abs(m5-n+k-4) + ...
            5*e^2*abs(m5-n+k-2)/2/(1-e^2).*b.^abs(m5-n+k-3) + ...
            5*e^3/2/(1-e^2)^(3/2).*b.^abs(m5-n+k-2);
        
        dg2de = abs(m2-n+k-2).*b.^(abs(m2-n+k-2))/e/sqrt(1-e^2);
        dg3de = abs(m3-n+k-2).*(abs(m3-n+k-3).*b.^(abs(m3-n+k-3)) + ...
            e/sqrt(1-e^2)*b.^(abs(m3-n+k-2)))/e/sqrt(1-e^2) + ...
            b.^abs(m3-n+k-2)/(1-e^2)^(3/2);
        dg4de = abs(m4-n+k-2).*(3*e^2/2/(1-e^2)*b.^(abs(m4-n+k-2)) +...
            abs(m4-n+k-3).*(3*e/2/sqrt(1-e^2)*b.^(abs(m4-n+k-3)) +...
            abs(m4-n+k-4)/2.*b.^abs(m4-n+k-4)))/e/sqrt(1-e^2) +...
            3/2*abs(m4-n+k-2)/(1-e^2)^(3/2).*b.^abs(m4-n+k-3) +...
            3*e/(1-e^2)^2*b.^abs(m4-n+k-2);
        dg5de = abs(m5-n+k-2).*(5*e^3/2/(1-e^2)^(3/2)*b.^abs(m5-n+k-2) + ...
            abs(m5-n+k-3).*(5*e^2/2/(1-e^2)*b.^abs(m5-n+k-3) + ...
            abs(m5-n+k-4).*(e/sqrt(1-e^2)*b.^abs(m5-n+k-4) + ...
            abs(m5-n+k-5)/6.*b.^abs(m5-n+k-5))))/e/sqrt(1-e^2) + ...
            abs(m5-n+k-3).*abs(m5-n+k-2)/(1-e^2)^(3/2).*b.^abs(m5-n+k-4) +...
            5*abs(m5-n+k-2)*e/(1-e^2)^2.*b.^abs(m5-n+k-3) +...
            15*e^2/2/(1-e^2)^(5/2)*b.^abs(m5-n+k-2);
        
        Jn = besselj(n,k*e);
        dJnde = 0.5*(besselj(n-1,k*e) - besselj(n+1,k*e));
        
        dAkda = dAkda + Jn*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dAkde = dAkde + k*dJnde*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]+...
            Jn*(dCde.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] + ...
            C.'*[da2de.'*g2+a2.'*dg2de; da3de.'*g3+a3.'*dg3de; da4de.'*g4+a4.'*dg4de; da5de.'*g5+a5.'*dg5de]);
        dAkdi = dAkdi + Jn*dCdi.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dAkdo = dAkdo + Jn*dCdo.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
        dAkdl = dAkda;
        
        dBkda = dBkda + Jn*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        dBkde = dBkde + k*dJnde*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5]+...
            Jn*(dSde.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] + ...
            S.'*[db1de.'*g2+b1.'*dg2de; db2de.'*g3+b2.'*dg3de; db3de.'*g4+b3.'*dg4de; db4de.'*g5+b4.'*dg5de]);
        dBkdi = dBkdi + Jn*dSdi.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        dBkdo = dBkdo + Jn*dSdo.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
        dBkdl = dBkda;
        
        Akda = Akda + dAkda;
        Akde = Akde + dAkde;
        Akdi = Akdi + dAkdi;
        Akdo = Akdo + dAkdo;
        Akdl = Akdl + dAkdl;
        
        Bkda = Bkda + dBkda;
        Bkde = Bkde + dBkde;
        Bkdi = Bkdi + dBkdi;
        Bkdo = Bkdo + dBkdo;
        Bkdl = Bkdl + dBkdl;
        
        n = n+1;        
    end
    dJ2daFreq(:,k+1) = dRda*[Akda;Bkda];
    dJ2deFreq(:,k+1) = R*[Akde;Bkde];
    dJ2diFreq(:,k+1) = R*[Akdi;Bkdi];
    dJ2doFreq(:,k+1) = R*[Akdo;Bkdo];
    dJ2dlFreq(:,k+1) = R*[k*Bkdl;-k*Akdl];
    
    k = k+1;
end
k = 0:kMax;
trigMat = [cos(k*M);sin(k*M)];

dJ2da = sum(dJ2daFreq.*trigMat,'all');
dJ2de = sum(dJ2deFreq.*trigMat,'all');
dJ2di = sum(dJ2diFreq.*trigMat,'all');
dJ2do = sum(dJ2doFreq.*trigMat,'all');
dJ2dl = sum(dJ2dlFreq.*trigMat,'all');

dOe = zeros(6,1);
n = sqrt(primary.mu/a^3);
eta = sqrt(1-e^2);

dOe(1) = 2/n/a*dJ2dl;
dOe(2) = eta^2/n/a^2/e*dJ2dl - eta/n/a^2/e*dJ2do;
dOe(3) = cos(i)/n/a^2/eta/sin(i)*dJ2do;
dOe(4) = 1/n/a^2/eta/sin(i)*dJ2di;
dOe(5) = eta/n/a^2/e*dJ2de - cos(i)/n/a^2/eta/sin(i)*dJ2di;
dOe(6) = n -2/n/a*dJ2da -eta^2/n/a^2/e*dJ2de;

dJ2 = [dJ2da; dJ2de; dJ2di; dJ2do; dJ2dl];