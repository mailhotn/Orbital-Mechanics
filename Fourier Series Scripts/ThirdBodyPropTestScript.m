clear
%% Setup
ic = [40000, 0.01, 10, 40, 30, 10].';
primary = Earth;
third = {Moon};
Sat = SingleSat(ic,primary,third);
Prop = Propagator(Sat);
T = 2*pi*sqrt(ic(1)^3/primary.mu);
t = 0:100:T*10;

%% Propagate

% Prop using Good solution
tic
[~,X1] = Prop.PropEci3B(t,'Stable');
oe1 = eci2oe(X1.',[],primary,'me');
oe1(6,:) = unwrap(oe1(6,:)*pi/180)*180/pi;
stableTime = toc;
% Prop using Taylor solution
tic
[~,X2] = Prop.PropEci3B(t,'Taylor');
oe2 = eci2oe(X2.',[],primary,'me');
oe2(6,:) = unwrap(oe2(6,:)*pi/180)*180/pi;
taylorTime = toc;
% Prop using Naive solution
tic
[~,X3] = Prop.PropEci3B(t,'Direct');
oe3 = eci2oe(X3.',[],primary,'me');
oe3(6,:) = unwrap(oe3(6,:)*pi/180)*180/pi;
directTime = toc;
% Prop using Fourier solution
tic
[~,oeF] = Prop.PropOeFourier3B(t,6);
oeF = oeF.';
fourTime = toc;

%% Calculate relative errors
err2 = abs(oe1-oe2)./oe1;
err3 = abs(oe1-oe3)./oe1;
errF = abs(oeF-oe1)./oe1;

%% Plot
lWidth = 1.5;
t = t/86400;
figure(1)
tiledlayout(3,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(t,oe1(1,:),t,oe2(1,:),t,oe3(1,:),'--',t,oeF(1,:),LineWidth=lWidth)
ylabel('a [km]')
xlabel('t [day]')
nexttile
plot(t,oe1(2,:),t,oe2(2,:),t,oe3(2,:),'--',t,oeF(2,:),LineWidth=lWidth)
ylabel('e')
xlabel('t [day]')
nexttile
plot(t,oe1(3,:),t,oe2(3,:),t,oe3(3,:),'--',t,oeF(3,:),LineWidth=lWidth)
ylabel('i [deg]')
xlabel('t [day]')
nexttile
plot(t,oe1(4,:),t,oe2(4,:),t,oe3(4,:),'--',t,oeF(4,:),LineWidth=lWidth)
ylabel('\Omega [deg]')
xlabel('t [day]')
nexttile
plot(t,oe1(5,:),t,oe2(5,:),t,oe3(5,:),'--',t,oeF(5,:),LineWidth=lWidth)
ylabel('\omega [deg]')
xlabel('t [day]')
nexttile
plot(t,oe1(6,:),t,oe2(6,:),t,oe3(6,:),'--',t,oeF(6,:),LineWidth=lWidth)
ylabel('M [deg]')
xlabel('t [day]')

figure(2)
tiledlayout(3,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(t,err2(1,:),t,err3(1,:),'--',t,errF(1,:),LineWidth=lWidth)
ylabel('a')
xlabel('t [day]')
nexttile
plot(t,err2(2,:),t,err3(2,:),'--',t,errF(2,:),LineWidth=lWidth)
ylabel('e')
xlabel('t [day]')
nexttile
plot(t,err2(3,:),t,err3(3,:),'--',t,errF(3,:),LineWidth=lWidth)
ylabel('i')
xlabel('t [day]')
nexttile
plot(t,err2(4,:),t,err3(4,:),'--',t,errF(4,:),LineWidth=lWidth)
ylabel('\Omega')
xlabel('t [day]')
nexttile
plot(t,err2(5,:),t,err3(5,:),'--',t,errF(5,:),LineWidth=lWidth)
ylabel('\omega')
xlabel('t [day]')
nexttile
plot(t,err2(6,:),t,err3(6,:),'--',t,errF(6,:),LineWidth=lWidth)
ylabel('M')
xlabel('t [day]')
figure(3)
bar([stableTime,taylorTime,directTime,fourTime])
