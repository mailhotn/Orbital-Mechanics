clear
%% Setup
primary = Earth;
third = {Moon,Sun};
% T = 86164/2;
% sma = ((T/2/pi)^2*primary.mu)^(1/3);
% ic = [sma, 0.74, 63.4, 30, 270, 10].'; % Molniya

T = 86164;
sma = ((T/2/pi)^2*primary.mu)^(1/3);
ic = [sma, 0.001, 5, 30, 30, 10].'; % GEO

Sat = SingleSat(ic,primary,third);
Prop = Propagator(Sat);
t = 0:1000:86400*200;
kMax = 4;
%% Propagate

% Prop using Good solution
tic
[~,X1] = Prop.PropEci3B(t,'Stable');
oe1 = eci2oe(X1.',[],primary,'me');
oe1(6,:) = unwrap(oe1(6,:)*pi/180)*180/pi - sqrt(primary.mu./oe1(1,:).^3).*t*180/pi;
stableTime = toc;

lan1 = oe1(5,:) + oe1(6,:)+oe1(4,:);
oe1(6,:) = lan1;
% Prop using Taylor solution
tic
[~,X2] = Prop.PropEci3B(t,'Taylor');
oe2 = eci2oe(X2.',[],primary,'me');
oe2(6,:) = unwrap(oe2(6,:)*pi/180)*180/pi- sqrt(primary.mu./oe2(1,:).^3).*t*180/pi;
taylorTime = toc;
lan2 = oe2(5,:) + oe2(6,:)+oe2(4,:);
% oe2(6,:) = lan2;
% Prop using Naive solution
tic
[~,X3] = Prop.PropEci3B(t,'Direct');
oe3 = eci2oe(X3.',[],primary,'me');
oe3(6,:) = unwrap(oe3(6,:)*pi/180)*180/pi- sqrt(primary.mu./oe3(1,:).^3).*t*180/pi;
directTime = toc;
lan3 = oe3(5,:) + oe3(6,:)+oe3(4,:);
% oe3(6,:) = lan3;
% Prop using Fourier solution
tic
[~,oeF] = Prop.PropOeFourier3B(t,kMax);
oeF = oeF.';
fourTime = toc;
% oeF(6,:) = me2ta(oeF(6,:),oeF(2,:));
oeF(6,:) = oeF(6,:) - sqrt(primary.mu./oeF(1,:).^3).*t*180/pi;
lanF = oeF(5,:) + oeF(6,:)+ oeF(4,:);
oeF(6,:) = lanF;
% Prop using singly averaged
tic
[~,oeS] = Prop.PropOeMeanSingle3B(t);
oeS = oeS.';
singleTime = toc;
oeS(6,:) = oeS(6,:) - sqrt(primary.mu./oeS(1,:).^3).*t*180/pi;
lanS =  oeS(5,:) + oeS(6,:) + oeS(4,:);
oe2 = oeS;
lan2 = lanS;
oe2(6,:) = lan2;
%% Calculate relative errors
err2 = abs(oe2-oe1)./abs(oe1);
err3 = abs(oe3-oe1)*0;
errF = abs(oeF-oe1)./abs(oe1);
meanErr = trapz(t,errF,2)/t(end)
%% Plot
lWidth = 1.5;
tD = t/86400;
figure(1)
tiledlayout(3,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(tD,oe1(1,:),tD,oe2(1,:),tD,oeF(1,:),'--',LineWidth=lWidth)
ylabel('a [km]')
xlabel('t [day]')

nexttile
plot(tD,oe1(2,:),tD,oe2(2,:),tD,oeF(2,:),'--',LineWidth=lWidth)
ylabel('e')
xlabel('t [day]')
legend('Numerical - Stable','Numerical - Taylor','Fourier')
nexttile
plot(tD,oe1(3,:),tD,oe2(3,:),tD,oeF(3,:),'--',LineWidth=lWidth)
ylabel('i [deg]')
xlabel('t [day]')
nexttile
plot(tD,oe1(4,:),tD,oe2(4,:),tD,oeF(4,:),'--',LineWidth=lWidth)
ylabel('\Omega [deg]')
xlabel('t [day]')
nexttile
plot(tD,oe1(5,:),tD,oe2(5,:),tD,oeF(5,:),'--',LineWidth=lWidth)
ylabel('\omega [deg]')
xlabel('t [day]')
nexttile
plot(tD,oe1(6,:),tD,oe2(6,:),tD,oeF(6,:),'--',LineWidth=lWidth)
ylabel('M - nt [deg]')
xlabel('t [day]')

figure(2)
tiledlayout(3,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(tD,err3(1,:),tD,err2(1,:),'--',tD,errF(1,:),LineWidth=lWidth)
ylabel('a')
xlabel('t [day]')
nexttile
plot(tD,err3(2,:),tD,err2(2,:),'--',tD,errF(2,:),LineWidth=lWidth)
ylabel('e')
xlabel('t [day]')
nexttile
plot(tD,err3(3,:),tD,err2(3,:),'--',tD,errF(3,:),LineWidth=lWidth)
ylabel('i')
xlabel('t [day]')
nexttile
plot(tD,err3(4,:),tD,err2(4,:),'--',tD,errF(4,:),LineWidth=lWidth)
ylabel('\Omega')
xlabel('t [day]')
nexttile
plot(tD,err3(5,:),tD,err2(5,:),'--',tD,errF(5,:),LineWidth=lWidth)
ylabel('\omega')
xlabel('t [day]')
nexttile
plot(tD,err3(6,:),tD,err2(6,:),'--',tD,errF(6,:),LineWidth=lWidth)
ylabel('M - nt')
xlabel('t [day]')
figure(3)
bar([stableTime,taylorTime,directTime,fourTime,singleTime])

figure(4)
plot(tD,lan1,tD,lan2,tD,lanF,'--',LineWidth=lWidth)
ylabel('L - nt [deg]')
xlabel('t [day]')
legend('Numerical - Stable','Numerical - Taylor','Fourier')
