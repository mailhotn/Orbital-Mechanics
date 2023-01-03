clear
oe = [10000, 0.05, 40, 40, 10, 0];

Sat = SingleSat(oe,earth());
Prop = Propagator(Sat);
T = 2*pi*sqrt(oe(1)^3/Sat.primary.mu);
nOrb = 200;
t = 0:60:nOrb*T;
kLow = 4;
kHigh = 4;
%% Conventional ECI
[~,Xeci] = Prop.PropEciJ2(t);
oeC = eci2oe(Xeci(:,1:3),Xeci(:,4:6));
oeC(6,:) = ta2me(oeC(6,:),oeC(2,:));
%% Conventional LPE
tic
[~,oeC] = Prop.PropOeOsc(t);
oeC = oeC.';
oscTime = toc
%% Mean
tic
[~,OeM] = Prop.PropOeMeanFast(t);
oeB = me2osc(OeM.');
bTime = toc
%% Mean SP only
[~,OeS] = Prop.PropOeMeanShort(t);
oeS = me2oscSP(OeS.');
%% Fourier
% [~,oeF2] = Prop.PropOeFourier2(t,2);
% [~,oeF5] = Prop.PropOeFourier2(t,5);
% tic
% [~,oeFLow] = Prop.PropOeFourier2(t,kLow);
% oeFLow = oeFLow.';
% F10Time = toc

tic
[~,oeFHigh] = Prop.PropOeFourier(t,kHigh);
oeFHigh = oeFHigh.';
F15Time = toc
t = t/T;
%% Plot
% oeFHigh = oeS;
% semimajor axis
f = figure(1);
% f.Position(3:4) = [500,250];
plot(t,oeC(1,:)/oe(1),t,oeFHigh(1,:)/oe(1),'--','linewidth',1.5)
legend('Conventional','Fourier')
ylabel('$\frac{a_x\left(t\right)}{a\left(0\right)}$','fontsize',18,'interpreter','latex')
xlabel('$Orbit$','interpreter','latex','fontsize',16)
grid on
xlim([0,10])
xticks(1:nOrb)

% % eccentricity
figure(2)
% f.Position(3:4) = [500,250];
plot(t,oeC(2,:),t,oeFHigh(2,:),'--','linewidth',1.5)
legend('Conventional','Fourier')
ylabel('$e_x$','fontsize',16,'interpreter','latex')
xlabel('$Orbit$','interpreter','latex','fontsize',16)
grid on
xlim([0,10])
xticks(1:nOrb)

% inclination
f = figure(3);
% f.Position(3:4) = [500,250];
plot(t,oeC(3,:),t,oeFHigh(3,:),'--','linewidth',1.5)
legend('Conventional','Fourier')
ylabel('$i_x \left[{rad}\right]$','fontsize',16,'interpreter','latex')
xlabel('$Orbit$','interpreter','latex','fontsize',16)
grid on
xlim([0,nOrb])
xticks(1:nOrb)

% % RAAN
figure(4)
% f.Position(3:4) = [500,250];
plot(t,oeC(4,:),t,oeFHigh(4,:),'--','linewidth',1.5)
legend('Conventional','Fourier')
ylabel('$\Omega_x \left[{rad}\right]$','fontsize',16,'interpreter','latex')
xlabel('$Orbit$','interpreter','latex','fontsize',16)
grid on
xlim([0,nOrb])
xticks(1:nOrb)

% % AOP
figure(5)
% f.Position(3:4) = [500,250];
plot(t,oeC(5,:),t,oeFHigh(5,:),'--','linewidth',1.5)
legend('Conventional','Fourier')
ylabel('$\omega_x \left[{rad}\right]$','fontsize',16,'interpreter','latex')
xlabel('$Orbit$','interpreter','latex','fontsize',16)
grid on
xlim([0,nOrb])
xticks(1:nOrb)

% % mean anomaly
figure(6)
% f.Position(3:4) = [500,250];
plot(t,oeC(6,:),t,oeFHigh(6,:),'--','linewidth',1.5)
legend('Conventional','Fourier')
ylabel('$M_x \left[{rad}\right]$','fontsize',16,'interpreter','latex')
xlabel('$Orbit$','interpreter','latex','fontsize',16)
grid on
xlim([0,nOrb])
xticks(1:nOrb)



%% Plot errors
errB = abs(oeC-oeB);
% errFLow = abs(oeC-oeFLow);
errFHigh = abs(oeC-oeFHigh);
errB = [errB(1,:)/oe(1);errB(2,:)/oe(2);errB(3:end,:)*pi/180];
% errFLow = [errFLow(1,:)/oe(1);errFLow(2,:)/oe(2);errFLow(3:end,:)*pi/180];
errFHigh = [errFHigh(1,:)/oe(1);errFHigh(2,:)/oe(2);errFHigh(3:end,:)*pi/180];
% intErrB = trapz(t.',errB,2)/t(end)
% intErr10 = trapz(t.',errFLow,2)/t(end)
intErr15 = trapz(t.',errFHigh,2)/t(end)

% labelHigh = ['Fourier'];
% 
% % 
% figure(11) % semimajor axis
% plot(t,errFHigh(1,:),t,errB(1,:),'linewidth',2)
% legend(labelHigh,'Brouwer','fontsize',12,'location','best')
% ylabel('$\frac{\left|a_c-a_x\right|}{a\left(0\right)}$','fontsize',18,'interpreter','latex')
% xlabel('$Orbit$','interpreter','latex','fontsize',18)
% grid on
% xlim([0,nOrb])
% xticks(1:nOrb)
% figure(12) % eccentricity
% plot(t,errFHigh(2,:),t,errB(2,:),'linewidth',2)
% legend(labelHigh,'Brouwer','fontsize',12,'location','best')
% ylabel('$\frac{\left|e_c-e_x\right|}{e\left(0\right)}$','fontsize',18,'interpreter','latex')
% xlabel('$Orbit$','interpreter','latex','fontsize',18)
% grid on
% xlim([0,nOrb])
% xticks(1:nOrb)
% figure(13) % inclination
% plot(t,errFHigh(3,:),t,errB(3,:),'linewidth',2)
% legend(labelHigh,'Brouwer','fontsize',12,'location','best')
% ylabel('$\left|i_c-i_x\right| \left[{rad}\right]$','fontsize',18,'interpreter','latex')
% xlabel('$Orbit$','interpreter','latex','fontsize',18)
% grid on
% xlim([0,nOrb])
% xticks(1:nOrb)
% figure(14) % raan
% plot(t,errFHigh(4,:),t,errB(4,:),'linewidth',2)
% legend(labelHigh,'Brouwer','fontsize',12,'location','best')
% ylabel('$\left|\Omega_c-\Omega_x\right| \left[{rad}\right]$','fontsize',18,'interpreter','latex')
% xlabel('$Orbit$','interpreter','latex','fontsize',18)
% grid on
% xlim([0,nOrb])
% xticks(1:nOrb)
% figure(15) % aop
% plot(t,errFHigh(5,:),t,errB(5,:),'linewidth',2)
% legend(labelHigh,'Brouwer','fontsize',12,'location','best')
% ylabel('$\left|\omega_c-\omega_x\right| \left[{rad}\right]$','fontsize',18,'interpreter','latex')
% xlabel('$Orbit$','interpreter','latex','fontsize',18)
% grid on
% xlim([0,nOrb])
% xticks(1:nOrb)
% figure(16) % mean anomaly
% plot(t,errFHigh(6,:),t,errB(6,:),'linewidth',2)
% legend(labelHigh,'Brouwer','fontsize',12,'location','best')
% ylabel('$\left|M_c-M_x\right| \left[{rad}\right]$','fontsize',18,'interpreter','latex')
% xlabel('$Orbit$','interpreter','latex','fontsize',18)
% grid on
% xlim([0,nOrb])
% if (max(errFLow(6,:))>0.2) || (max(errB(6,:))>0.2)
%     ylim([0,0.2])
% end
% xticks(1:nOrb)
