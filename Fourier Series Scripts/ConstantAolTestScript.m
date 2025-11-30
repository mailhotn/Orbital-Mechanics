primary = earth();
J2 = primary.J2;
Re = primary.Re;
mu = primary.mu;

a = 10000;
eRange = linspace(0,0.3,1000);
iRange = linspace(0,90,1000);

[I,E] = meshgrid(iRange,eRange);
eta = sqrt(1-E.^2);
J1 = besselj(1,E);
dJ1de = 0.5*(besselj(0,E) -besselj(2,E));

% Average element rate due to J2
dAop = 0.75*J2*sqrt(mu)*Re^2*a^(-7/2)*(1-E.^2).^-2.*(5*cosd(I).^2-1);
dMan = 0.75*J2*sqrt(mu)*Re^2*a^(-7/2)*(1-E.^2).^(-3/2).*(3*cosd(I).^2-1);
dRan = -1.5*J2*sqrt(mu)*Re^2*a^(-7/2)*(1-E.^2).^-2.*cosd(I);

% Average element rate due to ur = cos(M) i.e. 1 km/s thrust magnitude
% dAolCont = 0.5*sqrt(a/mu)/e*(-2*eta^3/e*J1 +(2*eta^4/e*J1+4*e^4*dJ1de))
dAopCont = 0.5*sqrt(a/mu)*(-2*eta.^3./E.^2.*J1);

constAopRanCont = abs(dAop+dRan)./abs(dAopCont)*1000; % thrust m/s^2 for constant aop + raan
thrustCutoff = 1e-4; % m/s
constAopRanCont(constAopRanCont>thrustCutoff) = thrustCutoff; % cutoff
% figure(1)
% plot(iRange,(dAop+dMan)./abs(dAopCont))
figure(1)
contourf(I,E,abs(dAop+dRan))
colorbar


figure(2)
shading interp
contourf(I,E,constAopRanCont*86400)
colorbar
colormap jet
hold on
contour(I,E,dAop+dRan,[0,0],'w')
xlabel('i [deg]')
ylabel('e')
hold off

