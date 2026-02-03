function dX = ConDynOeOsc(P,t,X)
% Controlled, GVE equations & J2 forces
OE  = reshape(X,6,P.Con.nSats);
% Normalized Parameters
mu = 1;
Re = 1;
% Element vectors (angles already in radians)
sma   = OE(1,:);
ecc   = OE(2,:);
inc = OE(3,:);
aop   = OE(5,:);
man  = OE(6,:);
tan  = pi/180*me2ta(man*180/pi,ecc); %  conversion M to f
% Secondary definitions
p = sma.*(1-ecc.^2);
h = sqrt(mu*p);
r = p./(1+ecc.*cos(tan));
n = sqrt(mu./sma.^3);
aol = tan + aop;
% GVE
MJ2 = -3*mu./r.^4.*P.Con.J2.*Re.^2; % common coefficient
% Forces - J2 + Controller
fCont = P.Control.ControlRSW(t,OE);
fR = MJ2/2.*(1-3*sin(inc).^2.*sin(aol).^2) + fCont(1,:);
fS = MJ2.*sin(inc).^2.*sin(aol).*cos(aol) + fCont(2,:);
fH = MJ2.*cos(inc).*sin(inc).*sin(aol) + fCont(3,:); 
% Element Rates
da = 2*sma.^2./h.*ecc.*sin(tan).*fR +...
    2*sma.^2./h.*(1 + ecc.*cos(tan)).*fS;
de = p./h.*sin(tan).*fR +...
    r./h.*(ecc + 2*cos(tan) + ecc.*cos(tan).^2).*fS;
di = r./h.*cos(aol).*fH;
dO = r.*sin(aol)./(h.*sin(inc)).*fH;
dw = -p./(h.*ecc).*cos(tan).*fR +...
    r./(h.*ecc).*(2+ecc.*cos(tan)).*sin(tan).*fS +...
    -r./h.*sin(aol).*cos(inc)./sin(inc).*fH;
dM = (p.*cos(tan)-2*r.*ecc)./(n.*sma.^2.*ecc).*fR +...
    -(p+r).*sin(tan)./(n.*sma.^2.*ecc).*fS;

dOe = reshape([da;de;di;dO;dw;dM + n],6*P.Con.nSats,1);
% Equations of Motion
dX = dOe;
end