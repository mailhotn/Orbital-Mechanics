clear
oe = [10000, 0.1, 63.25, 30, 70, 50];

Sat = SingleSat(oe,earth());
Prop = Propagator(Sat);
T = 2*pi*sqrt(oe(1)^3/Sat.primary.mu);
nOrb = 2;
t = 0:60:nOrb*T;
kLow = 4;
kHigh = 5;
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
%% Fourier
[~,oeF2] = Prop.PropOeFourier2(t,2);
[~,oeF5] = Prop.PropOeFourier2(t,5);
tic
[~,oeFLow] = Prop.PropOeFourier2(t,kLow);
F10Time = toc
tic
[~,oeFHigh] = Prop.PropOeFourier2(t,kHigh);
F15Time = toc
t = t/T;
%% Plot
% % semimajor axis
% figure(1)
% plot(t,oeC(1,:),t,oeB(1,:),'.'...
%     ,t,oeF2(1,:),'--',t,oeF5(1,:),'--',t,oeFLow(1,:),'--',t,oeFHigh(1,:),'--','linewidth',1.5)
% legend('ECI','Brouwer','F2','F5','F10','F15')
% % eccentricity
% figure(2)
% plot(t,oeC(2,:),t,oeB(2,:),'.'...
%     ,t,oeF2(2,:),'--',t,oeF5(2,:),'--',t,oeFLow(2,:),'--',t,oeFHigh(2,:),'--','linewidth',1.5)
% legend('ECI','Brouwer','F2','F5','F10','F15')
% % inclination
% figure(3)
% plot(t,oeC(3,:),t,oeB(3,:),'.'...
%     ,t,oeF2(3,:),'--',t,oeF5(3,:),'--',t,oeFLow(3,:),'--',t,oeFHigh(3,:),'--','linewidth',1.5)
% legend('ECI','Brouwer','F2','F5','F10','F15')
% % RAAN
% figure(4)
% plot(t,oeC(4,:),t,oeB(4,:),'.'...
%     ,t,oeF2(4,:),'--',t,oeF5(4,:),'--',t,oeFLow(4,:),'--',t,oeFHigh(4,:),'--','linewidth',1.5)
% legend('ECI','Brouwer','F2','F5','F10','F15')
% % AOP
% figure(5)
% plot(t,oeC(5,:),t,oeB(5,:),'.'...
%     ,t,oeF2(5,:),'--',t,oeF5(5,:),'--',t,oeFLow(5,:),'--',t,oeFHigh(5,:),'--','linewidth',1.5)
% legend('ECI','Brouwer','F2','F5','F10','F15')
% % mean anomaly
% figure(6)
% plot(t,oeC(6,:),t,oeB(6,:),'.'...
%     ,t,oeF2(6,:),'--',t,oeF5(6,:),'--',t,oeFLow(6,:),'--',t,oeFHigh(6,:),'--','linewidth',1.5)
% legend('ECI','Brouwer','F2','F5','F10','F15')

%% Plot errors
errB = abs(oeC-oeB);
errFLow = abs(oeC-oeFLow);
errFHigh = abs(oeC-oeFHigh);
errB = [errB(1,:)/oe(1);errB(2,:)/oe(2);errB(3:end,:)*pi/180];
errFLow = [errFLow(1,:)/oe(1);errFLow(2,:)/oe(2);errFLow(3:end,:)*pi/180];
errFHigh = [errFHigh(1,:)/oe(1);errFHigh(2,:)/oe(2);errFHigh(3:end,:)*pi/180];
intErrB = trapz(t.',errB,2)/t(end)
intErr10 = trapz(t.',errFLow,2)/t(end)
intErr15 = trapz(t.',errFHigh,2)/t(end)

labelHigh = ['Fourier'];

% 
% figure(11) % semimajor axis
% plot(t,errB(1,:),t,errFHigh(1,:),'--','linewidth',2)
% legend('Brouwer',labelHigh,'fontsize',18)
% ylabel('$\frac{\left|a_c-a_x\right|}{a\left(0\right)}$','fontsize',18,'interpreter','latex')
% xlabel('$Orbit$','interpreter','latex','fontsize',18)
% grid on
% xlim([0,nOrb])
% xticks(1:nOrb)
figure(12) % eccentricity
plot(t,errFHigh(2,:),t,errB(2,:),'--','linewidth',2)
legend(labelHigh,'Brouwer','fontsize',12)
ylabel('$\frac{\left|e_c-e_x\right|}{e\left(0\right)}$','fontsize',18,'interpreter','latex')
xlabel('$Orbit$','interpreter','latex','fontsize',18)
grid on
xlim([0,nOrb])
xticks(1:nOrb)
% figure(13) % inclination
% plot(t,errB(3,:),t,errFHigh(3,:),'--','linewidth',2)
% legend('Brouwer',labelHigh)
% ylabel('$\left|i_c-i_x\right| \left[{rad}\right]$','fontsize',18,'interpreter','latex')
% xlabel('$Orbit$','interpreter','latex','fontsize',18)
% grid on
% xlim([0,nOrb])
% xticks(1:nOrb)
% figure(14) % raan
% plot(t,errB(4,:),t,errFHigh(4,:),'--','linewidth',2)
% legend('Brouwer',labelHigh)
% ylabel('$\left|\Omega_c-\Omega_x\right| \left[{rad}\right]$','fontsize',18,'interpreter','latex')
% xlabel('$Orbit$','interpreter','latex','fontsize',18)
% grid on
% xlim([0,nOrb])
% xticks(1:nOrb)
% figure(15) % aop
% plot(t,errB(5,:),t,errFHigh(5,:),'--','linewidth',2)
% legend('Brouwer',labelHigh)
% ylabel('$\left|\omega_c-\omega_x\right| \left[{rad}\right]$','fontsize',18,'interpreter','latex')
% xlabel('$Orbit$','interpreter','latex','fontsize',18)
% grid on
% xlim([0,nOrb])
% xticks(1:nOrb)
% figure(16) % mean anomaly
% plot(t(2:end),errB(6,2:end),t(2:end),errFHigh(6,2:end),'--','linewidth',2)
% legend('Brouwer',labelHigh)
% ylabel('$\left|M_c-M_x\right| \left[{rad}\right]$','fontsize',18,'interpreter','latex')
% xlabel('$Orbit$','interpreter','latex','fontsize',18)
% grid on
% xlim([0,nOrb])
% if (max(errFLow(6,:))>0.2) || (max(errB(6,:))>0.2)
%     ylim([0,0.2])
% end
% xticks(1:nOrb)
