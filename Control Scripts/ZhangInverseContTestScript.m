clear
%% Setup
icM = [7000,0.005,0.175*180/pi,10,10,0].';
tOe = [7191.1,0.005,0.175*180/pi].';
icOsc = me2osc(icM);

primary = Earth;
dScale = primary.Re;
tScale = sqrt(primary.Re^3/primary.mu);
fScale = dScale/tScale^2;

T = 2*pi*sqrt(icM(1)^3/primary.mu);
thrustMag = 6e-6;
tol = 1e-5;

Sat = SingleSat(icOsc,primary);
Control = ZhangInverse(tOe,thrustMag,tol);
Prop = Propagator(Sat,Control);

t = 0:10:T*10;

[~,X] = Prop.PropConOeOsc(t);
X = X.';
XM = osc2me(X);

%% Plot
lWidth = 1.5;
tD = t/T;
% X = XM;
figure(1)
tiledlayout(3,2,"Padding","compact","Tilespacing","tight")

nexttile
plot(tD,X(1,:),tD,XM(1,:),'--',LineWidth=lWidth)
ylabel('a [km]')
xlabel('t [orbit]')
nexttile
plot(tD,X(2,:),tD,XM(2,:),'--',LineWidth=lWidth)
ylabel('e')
xlabel('t [orbit]')
legend('Osculating','Mean')
nexttile
plot(tD,X(3,:),tD,XM(3,:),'--',LineWidth=lWidth)
ylabel('i [deg]')
xlabel('t [orbit]')
nexttile
plot(tD,X(4,:),tD,XM(4,:),'--',LineWidth=lWidth)
ylabel('\Omega [deg]')
xlabel('t [orbit]')
nexttile
plot(tD,X(5,:),tD,XM(5,:),'--',LineWidth=lWidth)
ylabel('\omega [deg]')
xlabel('t [orbit]')
nexttile
plot(tD,X(6,:),tD,XM(6,:),'--',LineWidth=lWidth)
ylabel('M [deg]')
xlabel('t [orbit]')


f = Control.ControlRSW([],X);

figure(2)
plot(t/T,f(1,:))