clear
oe = [10000, 0.3, 50, 30, 70, 0];

Sat = SingleSat(oe,earth());
Prop = Propagator(Sat);
t = 0:100:86400;
%% Conventional
[~,Xeci] = Prop.PropEciJ2(t);
oeC = eci2oe(Xeci(:,1:3),Xeci(:,4:6));
oeC(6,:) = ta2me(oeC(6,:),oeC(2,:));
%% Mean
[~,OeM] = Prop.PropOeMeanFast(t);
oeB = me2osc(OeM.');
%% Fourier
[~,oeF2] = Prop.PropOeFourier(t,2);
[~,oeF5] = Prop.PropOeFourier(t,5);
[~,oeF10] = Prop.PropOeFourier(t,10);
[~,oeF15] = Prop.PropOeFourier(t,15);

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
