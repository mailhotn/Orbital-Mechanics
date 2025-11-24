primary = earth();
J2 = primary.J2;
Re = primary.Re;
mu = primary.mu;

a = 7000;
e = 0.01;
eta = sqrt(1-e^2);
iRange = linspace(0,90,1000);


J1 = besselj(1,e);
dJ1de = 0.5*(besselj(0,e) -besselj(2,e));

% Average element rate due to J2
dAop = 0.75*J2*sqrt(mu)*Re^2*a^(-7/2)*(1-e^2)^-2*(5*cosd(iRange).^2-1);
dMan = 0.75*J2*sqrt(mu)*Re^2*a^(-7/2)*(1-e^2)^(-3/2)*(3*cosd(iRange).^2-1);

% Average element rate due to ur = cos(M) i.e. 1 km/s thrust magnitude
dAolCont = 0.5*sqrt(a/mu)/e*(-2*eta^3/e*J1 +(2*eta^4/e*J1+4*e^4*dJ1de))

figure(1)
plot(iRange,(dAop+dMan)./abs(dAolCont))
