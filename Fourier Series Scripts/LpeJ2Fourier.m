function dOe = LpeJ2Fourier(t,oe,kMax,primary)

if nargin < 3
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
dR0di = 3*cos(i)*sin(i)/a^3/(1-e^2)^(3/2);
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


C1xA1 = 6*(3*e*sin(i)^2*cos(aop)^2-1)/(1-e^2)^2;

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

dSdi = [12*sin(i)*cos(i)*sin(aop)*cos(aop)/2/(1-e^2)^3;
    36*e*sin(i)*cos(i)*sin(aop)*cos(aop)/4/(1-e^2)^(7/2);
    36*e^2*sin(i)*cos(i)*sin(aop)*cos(aop)/8/(1-e^2)^4;
    12*e^3*sin(i)*cos(i)*sin(aop)*cos(aop)/16/(1-e^2)^(9/2)];

dSdo = [6*sin(i)^2*(2*cos(aop)^2-1)/2/(1-e^2)^3;
    18*e*sin(i)^2*(2*cos(aop)^2-1)/4/(1-e^2)^(7/2);
    18*e^2*sin(i)^2*(2*cos(aop)^2-1)/8/(1-e^2)^4;
    6*e^3*sin(i)^2*(2*cos(aop)^2-1)/16/(1-e^2)^(9/2)];

dSde = [-18*e*sin(i)^2*sin(aop)*cos(aop)/(1-e^2)^4;
         -9*sin(i)^2*sin(aop)*cos(aop)*(6*e^2+1)/2/(1-e^2)^(9/2);
         -9*e*sin(i)^2*sin(aop)*cos(aop)*(3*e^2+1)/2/(1-e^2)^5;
         -9*e^2*sin(i)^2*sin(aop)*cos(aop)*(2*e^2+1)/8/(1-e^2)^(11/2)];

dJ2daFreq = dRda*[R0;0];
dJ2deFreq = R*[dR0de;0];
dJ2diFreq = R*[dR0di;0];
dJ2doFreq = R*[dR0do;0];
dJ2dlFreq = R*[dR0dl;0];

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
    
    Akda = C1xA1*besselj(k,k*e) +...
        besselj(n,-k*e)*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
    Akde = (C1xA1*k*dJkde + dC1xA1de*Jk) + ...
        -k*dJnde*C.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5]+...
        Jn*(dCde.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5] + ...
        C.'*[da2de.'*g2+a2.'*dg2de; da3de.'*g3+a3.'*dg3de; da4de.'*g4+a4.'*dg4de; da5de.'*g5+a5.'*dg5de]);
    Akdi = dC1xA1di*besselj(k,k*e) +...
        besselj(n,-k*e)*dCdi.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
    Akdo = dC1xA1do*besselj(k,k*e) +...
        besselj(n,-k*e)*dCdo.'*[a2.'*g2; a3.'*g3; a4.'*g4; a5.'*g5];
    Akdl = Akda;
    
    Bkda = besselj(n,-k*e)*S.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
    Bkde = -k*dJnde*SS.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5]+...
        Jn*(dSde.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5] + ...
        S.'*[db1de.'*g2+b1.'*dg2de; db2de.'*g3+b2.'*dg3de; db3de.'*g4+b3.'*dg4de; db4de.'*g5+b4.'*dg5de]);
    Bkdi = besselj(n,-k*e)*dSdi.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
    Bkdo = besselj(n,-k*e)*dSdo.'*[b1.'*g2; b2.'*g3; b3.'*g4; b4.'*g5];
    Bkdl = Bkda;
    while n <= k + 5
        
    end
    k = k+1;
end
    
    