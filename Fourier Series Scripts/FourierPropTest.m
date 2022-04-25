clear
oe = [10000, 0.2, 63.0, 30, 70, 0];

Sat = SingleSat(oe,earth());
Prop = Propagator(Sat);
t = 0:100:86400;
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
[~,oeF10] = Prop.PropOeFourier2(t,10);
F10Time = toc
tic
[~,oeF15] = Prop.PropOeFourier2(t,15);
F15Time = toc

%% Plot
% semimajor axis
figure(1)
plot(t,oeC(1,:),t,oeB(1,:),'.'...
    ,t,oeF2(1,:),'--',t,oeF5(1,:),'--',t,oeF10(1,:),'--',t,oeF15(1,:),'--','linewidth',1.5)
legend('ECI','Brouwer','F2','F5','F10','F15')
% eccentricity
figure(2)
plot(t,oeC(2,:),t,oeB(2,:),'.'...
    ,t,oeF2(2,:),'--',t,oeF5(2,:),'--',t,oeF10(2,:),'--',t,oeF15(2,:),'--','linewidth',1.5)
legend('ECI','Brouwer','F2','F5','F10','F15')
% inclination
figure(3)
plot(t,oeC(3,:),t,oeB(3,:),'.'...
    ,t,oeF2(3,:),'--',t,oeF5(3,:),'--',t,oeF10(3,:),'--',t,oeF15(3,:),'--','linewidth',1.5)
legend('ECI','Brouwer','F2','F5','F10','F15')
% RAAN
figure(4)
plot(t,oeC(4,:),t,oeB(4,:),'.'...
    ,t,oeF2(4,:),'--',t,oeF5(4,:),'--',t,oeF10(4,:),'--',t,oeF15(4,:),'--','linewidth',1.5)
legend('ECI','Brouwer','F2','F5','F10','F15')
% AOP
figure(5)
plot(t,oeC(5,:),t,oeB(5,:),'.'...
    ,t,oeF2(5,:),'--',t,oeF5(5,:),'--',t,oeF10(5,:),'--',t,oeF15(5,:),'--','linewidth',1.5)
legend('ECI','Brouwer','F2','F5','F10','F15')
% mean anomaly
figure(6)
plot(t,oeC(6,:),t,oeB(6,:),'.'...
    ,t,oeF2(6,:),'--',t,oeF5(6,:),'--',t,oeF10(6,:),'--',t,oeF15(6,:),'--','linewidth',1.5)
legend('ECI','Brouwer','F2','F5','F10','F15')

%% Plot errors
errB = abs(oeC-oeB);
errF10 = abs(oeC-oeF10);
errF15 = abs(oeC-oeF15);
errB = [errB(1,:)/oe(1);errB(2,:)/oe(2);errB(3:end,:)*pi/180];
errF10 = [errF10(1,:)/oe(1);errF10(2,:)/oe(2);errF10(3:end,:)*pi/180];
errF15 = [errF15(1,:)/oe(1);errF15(2,:)/oe(2);errF15(3:end,:)*pi/180];

figure(11) % semimajor axis
plot(t,errB(1,:),t,errF10(1,:),t,errF15(1,:))
legend('Brouwer','F10','F15')
xlabel('$Time \left[s\right]$','interpreter','latex','fontsize',12)
ylabel('$\frac{\left|a_c-a_x\right|}{a\left(0\right)}$','fontsize',15,'interpreter','latex')
grid on
xlim([0,86400])
figure(12) % eccentricity
plot(t,errB(2,:),t,errF10(2,:),t,errF15(2,:))
legend('Brouwer','F10','F15')
xlabel('$Time \left[s\right]$','interpreter','latex','fontsize',12)
ylabel('$\frac{\left|e_c-e_x\right|}{e\left(0\right)}$','fontsize',15,'interpreter','latex')
grid on
xlim([0,86400])
figure(13) % inclination
plot(t,errB(3,:),t,errF10(3,:),t,errF15(3,:))
legend('Brouwer','F10','F15')
xlabel('$Time \left[s\right]$','interpreter','latex','fontsize',12)
ylabel('$\left|i_c-i_x\right| \left[{rad}\right]$','fontsize',12,'interpreter','latex')
grid on
xlim([0,86400])
figure(14) % raan
plot(t,errB(4,:),t,errF10(4,:),t,errF15(4,:))
legend('Brouwer','F10','F15')
xlabel('$Time \left[s\right]$','interpreter','latex','fontsize',12)
ylabel('$\left|\Omega_c-\Omega_x\right| \left[{rad}\right]$','fontsize',12,'interpreter','latex')
grid on
xlim([0,86400])
figure(15) % aop
plot(t,errB(5,:),t,errF10(5,:),t,errF15(5,:))
legend('Brouwer','F10','F15')
xlabel('$Time \left[s\right]$','interpreter','latex','fontsize',12)
ylabel('$\left|\omega_c-\omega_x\right| \left[{rad}\right]$','fontsize',12,'interpreter','latex')
grid on
xlim([0,86400])
figure(16) % mean anomaly
plot(t(2:end),errB(6,2:end),t(2:end),errF10(6,2:end),t(2:end),errF15(6,2:end))
legend('Brouwer','F10','F15')
xlabel('$Time \left[s\right]$','interpreter','latex','fontsize',12)
ylabel('$\left|M_c-M_x\right| \left[{rad}\right]$','fontsize',12,'interpreter','latex')
grid on
xlim([0,86400])
