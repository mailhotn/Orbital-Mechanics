clear
%%
saveFolder = 'G:\My Drive\Doc Publications\2026 3 Body Fourier Paper\Figures';
% saveFolder = []; % comment to overwrite saved images

primary = Earth;
third = {Moon,Sun};
kMax = 4;
dT = 1800; % 
nDay = 30;
t = 0:dT:86400*nDay;
%% Setup Molniya
T = 86164/2;
sma = ((T/2/pi)^2*primary.mu)^(1/3);
ic = [sma, 0.74, 63.4, 30, 270, 10].'; % Molniya

MolSat = SingleSat(ic,primary,third);
MolProp = Propagator(MolSat);
%% Propagate Molniya
% Prop using Good solution
tic
[~,X1] = MolProp.PropEci3B(t,'Stable');
oeN = eci2oe(X1.',[],primary,'me');
oeN(6,:) = unwrap(oeN(6,:)*pi/180)*180/pi - sqrt(primary.mu./oeN(1,:).^3).*t*180/pi;
stableTime = toc;
% lan1 = oeN(5,:) + oeN(6,:)+oeN(4,:);
% oeN(6,:) = lan1;

% Prop using Taylor solution
% tic
% [~,X2] = Prop.PropEci3B(t,'Taylor');
% oe2 = eci2oe(X2.',[],primary,'me');
% oe2(6,:) = unwrap(oe2(6,:)*pi/180)*180/pi- sqrt(primary.mu./oe2(1,:).^3).*t*180/pi;
% taylorTime = toc;
% lan2 = oe2(5,:) + oe2(6,:)+oe2(4,:);
% % oe2(6,:) = lan2;
% % Prop using Naive solution
% tic
% [~,X3] = Prop.PropEci3B(t,'Direct');
% oe3 = eci2oe(X3.',[],primary,'me');
% oe3(6,:) = unwrap(oe3(6,:)*pi/180)*180/pi- sqrt(primary.mu./oe3(1,:).^3).*t*180/pi;
% directTime = toc;
% lan3 = oe3(5,:) + oe3(6,:)+oe3(4,:);
% oe3(6,:) = lan3;

% Prop using Fourier solution
tic
[~,oeF] = MolProp.PropOeFourier3B(t,kMax);
oeF = oeF.';
fourTime = toc;
xF = oe2eci(oeF,Earth,'me');
% oeF(6,:) = me2ta(oeF(6,:),oeF(2,:));
oeF(6,:) = oeF(6,:) - sqrt(primary.mu./oeF(1,:).^3).*t*180/pi;
% lanF = oeF(5,:) + oeF(6,:)+ oeF(4,:);
% oeF(6,:) = lanF;

% Prop using singly averaged
tic
[~,oeA] = MolProp.PropOeMeanSingle3B(t);
oeA = oeA.';
singleTime = toc;
xA = oe2eci(oeA,Earth,'me');
oeA(6,:) = oeA(6,:) - sqrt(primary.mu./oeA(1,:).^3).*t*180/pi;
%% Import Stela Data
oeS = importdata('G:\My Drive\Doc Data\Stela 3B Sims\MolSimIcrf.mat').';
xS = oe2eci(oeS,Earth,'me');
xS(:,end) = nan(6,1);
oeS(6,:) = unwrap(oeS(6,:)*pi/180)*180/pi - sqrt(primary.mu./oeS(1,:).^3).*t*180/pi;
oeS(6,end) = nan; % weird outlier

%% Calculate relative errors - Molniya
errA = abs(oeA-oeN);
errS = abs(oeS-oeN);
errS(6,end) = errS(6,end-1);
errF = abs(oeF-oeN);
dF = vecnorm((xF(1:3,:)-X1(:,1:3).'),2,1);
dA = vecnorm((xA(1:3,:)-X1(:,1:3).'),2,1);
dS = vecnorm((xS(1:3,:)-X1(:,1:3).'),2,1);
vF = vecnorm((xF(4:6,:)-X1(:,4:6).'),2,1);
vA = vecnorm((xA(4:6,:)-X1(:,4:6).'),2,1);
vS = vecnorm((xS(4:6,:)-X1(:,4:6).'),2,1);
meanErrFMol = trapz(t,errF,2)/t(end)
meanErrSMol = trapz(t,errA,2)/t(end)
meanErrStelaMol = trapz(t,errS,2)/t(end)
molTimes = [stableTime,fourTime,singleTime]
%% Plot Molniya
% plot elements
lWidth = 1.5;
tD = t/86400;
figure("Units","centimeters","Position",[0,0,16,16])
tiledlayout(3,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(tD,oeF(1,:),tD,oeA(1,:),tD,oeS(1,:),':',tD,oeN(1,:),'--',LineWidth=lWidth)
ylabel('$a \;[\rm{km}]$',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,oeF(2,:),tD,oeA(2,:),tD,oeS(2,:),':',tD,oeN(2,:),'--',LineWidth=lWidth)
ylabel('$e$',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
legend('Fourier','Averaged','STELA','Numerical',Location='northwest')
nexttile
plot(tD,oeF(3,:),tD,oeA(3,:),tD,oeS(3,:),':',tD,oeN(3,:),'--',LineWidth=lWidth)
ylabel('$i\; [\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,oeF(4,:),tD,oeA(4,:),tD,oeS(4,:),':',tD,oeN(4,:),'--',LineWidth=lWidth)
ylabel('$\Omega \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,oeF(5,:),tD,oeA(5,:),tD,oeS(5,:),':',tD,oeN(5,:),'--',LineWidth=lWidth)
ylabel('$\omega \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,oeF(6,:),tD,oeA(6,:),tD,oeS(6,:),':',tD,oeN(6,:),'--',LineWidth=lWidth)
ylabel('$M_0 \; [\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)

if ~isempty(saveFolder)
exportgraphics(gcf,[saveFolder '\MolniyaElements.eps'],...
        "ContentType","vector","Resolution",600)
end

% plot errors
figure("Units","centimeters","Position",[0,0,16,16])
tiledlayout(3,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(tD,errF(1,:),tD,errA(1,:),'--',tD,errS(1,:),':',LineWidth=lWidth)
ylabel('$\delta a \;[\rm{km}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,errF(2,:),tD,errA(2,:),'--',tD,errS(2,:),':',LineWidth=lWidth)
ylabel('$\delta e $',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
legend('Fourier','Averaged','STELA',Location='northwest')
nexttile
plot(tD,errF(3,:),tD,errA(3,:),'--',tD,errS(3,:),':',LineWidth=lWidth)
ylabel('$\delta i \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,errF(4,:),tD,errA(4,:),'--',tD,errS(4,:),':',LineWidth=lWidth)
ylabel('$\delta \Omega \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,errF(5,:),tD,errA(5,:),'--',tD,errS(5,:),':',LineWidth=lWidth)
ylabel('$\delta \omega \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,errF(6,:),tD,errA(6,:),'--',tD,errS(6,:),':',LineWidth=lWidth)
ylabel('$\delta M_0 \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)

if ~isempty(saveFolder)
exportgraphics(gcf,[saveFolder '\MolniyaErrors.eps'],...
        "ContentType","vector","Resolution",600)
end

% plot dist vel errors
figure("Units","centimeters","Position",[0,0,16,6])
tiledlayout(1,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(tD,dF,tD,dA,'--',tD,dS,':',LineWidth=lWidth)
ylabel('$\delta d \;[\rm{km}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,vF,tD,vA,'--',tD,vS,':',LineWidth=lWidth)
ylabel('$\delta v \;\left[\rm{km}\,\rm{s}^{-1}\right]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
legend('Fourier','Averaged','STELA',Location='northwest')

if ~isempty(saveFolder)
exportgraphics(gcf,[saveFolder '\MolniyaPosVelErr.eps'],...
        "ContentType","vector","Resolution",600)
end
%% Setup GEO
T = 86164;
sma = ((T/2/pi)^2*primary.mu)^(1/3);
ic = [sma, 0.001, 5, 30, 30, 10].'; % GEO

GeoSat = SingleSat(ic,primary,third);
GeoProp = Propagator(GeoSat);
%% Propagate Geo
% Prop using Good solution
tic
[~,X1] = GeoProp.PropEci3B(t,'Stable');
oeN = eci2oe(X1.',[],primary,'me');
oeN(6,:) = unwrap(oeN(6,:)*pi/180)*180/pi - sqrt(primary.mu./oeN(1,:).^3).*t*180/pi;
stableTime = toc;
lanN = oeN(5,:) + oeN(6,:)+oeN(4,:);
oeN(4,:) = lanN;

% Prop using Fourier solution
tic
[~,oeF] = GeoProp.PropOeFourier3B(t,kMax);
oeF = oeF.';
fourTime = toc;
xF = oe2eci(oeF,Earth,'me');
oeF(6,:) = oeF(6,:) - sqrt(primary.mu./oeF(1,:).^3).*t*180/pi;
lanF = oeF(5,:) + oeF(6,:)+ oeF(4,:);
oeF(4,:) = lanF;

% Prop using singly averaged
tic
[~,oeA] = GeoProp.PropOeMeanSingle3B(t);
oeA = oeA.';
singleTime = toc;
xA = oe2eci(oeA,Earth,'me');
oeA(6,:) = oeA(6,:) - sqrt(primary.mu./oeA(1,:).^3).*t*180/pi;
lanA =  oeA(5,:) + oeA(6,:) + oeA(4,:);
oeA(4,:) = lanA;
%% Import Stela Data
oeS = importdata('G:\My Drive\Doc Data\Stela 3B Sims\GeoSimIcrf.mat').';
xS = oe2eci(oeS,Earth,'me');
xS(:,end) = nan(6,1);
oeS(6,:) = unwrap(oeS(6,:)*pi/180)*180/pi - sqrt(primary.mu./oeS(1,:).^3).*t*180/pi;
lanS = oeS(5,:) + oeS(6,:)+oeS(4,:);
lanS(end) = nan; % weird outlier
oeS(4,:) = lanS;
%% Calculate relative errors - GEO
errA = abs(oeA-oeN);
errS = abs(oeS-oeN);
errS(4,end) = errS(4,end-1);
errF = abs(oeF-oeN);
dF = vecnorm((xF(1:3,:)-X1(:,1:3).'),2,1);
dA = vecnorm((xA(1:3,:)-X1(:,1:3).'),2,1);
dS = vecnorm((xS(1:3,:)-X1(:,1:3).'),2,1);
vF = vecnorm((xF(4:6,:)-X1(:,4:6).'),2,1);
vA = vecnorm((xA(4:6,:)-X1(:,4:6).'),2,1);
vS = vecnorm((xS(4:6,:)-X1(:,4:6).'),2,1);
meanErrFGeo = trapz(t,errF(1:4,:),2)/t(end)
meanErrSGeo = trapz(t,errA(1:4,:),2)/t(end)
meanErrStelaGeo = trapz(t,errS(1:4,:),2)/t(end)
geoTimes = [stableTime,fourTime,singleTime]
%% Plot GEO
% plot elements
lWidth = 1.5;
tD = t/86400;
figure("Units","centimeters","Position",[0,0,16,12])
tiledlayout(2,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(tD,oeF(1,:),tD,oeA(1,:),tD,oeS(1,:),':',tD,oeN(1,:),'--',LineWidth=lWidth)
ylabel('$a \;[\rm{km}]$',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,oeF(2,:),tD,oeA(2,:),tD,oeS(2,:),':',tD,oeN(2,:),'--',LineWidth=lWidth)
ylabel('$e$',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,oeF(3,:),tD,oeA(3,:),tD,oeS(3,:),':',tD,oeN(3,:),'--',LineWidth=lWidth)
ylabel('$i\; [\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,lanF,tD,lanA,tD,lanS,':',tD,lanN,'--',LineWidth=lWidth)
ylabel('$L \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
legend('Fourier','Averaged','STELA','Numerical',Location='southwest')

if ~isempty(saveFolder)
exportgraphics(gcf,[saveFolder '\GeoElements.eps'],...
        "ContentType","vector","Resolution",600)
end

% plot errors
figure("Units","centimeters","Position",[0,0,16,12])
tiledlayout(2,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(tD,errF(1,:),tD,errA(1,:),'--',tD,errS(1,:),':',LineWidth=lWidth)
ylabel('$\delta a \;[\rm{km}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,errF(2,:),tD,errA(2,:),'--',tD,errS(2,:),':',LineWidth=lWidth)
ylabel('$\delta e $',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
legend('Fourier','Averaged','STELA',Location='northwest')
nexttile
plot(tD,errF(3,:),tD,errA(3,:),'--',tD,errS(3,:),':',LineWidth=lWidth)
ylabel('$\delta i \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,errF(4,:),tD,errA(4,:),'--',tD,errS(4,:),':',LineWidth=lWidth)
ylabel('$\delta L \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)

if ~isempty(saveFolder)
exportgraphics(gcf,[saveFolder '\GeoErrors.eps'],...
        "ContentType","vector","Resolution",600)
end

% plot dist vel errors
figure("Units","centimeters","Position",[0,0,16,6])
tiledlayout(1,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(tD,dF,tD,dA,'--',tD,dS,':',LineWidth=lWidth)
ylabel('$\delta d \;[\rm{km}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,vF,tD,vA,'--',tD,vS,':',LineWidth=lWidth)
ylabel('$\delta v \;\left[\rm{km}\,\rm{s}^{-1}\right]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
legend('Fourier','Averaged','STELA',Location='northwest')

if ~isempty(saveFolder)
exportgraphics(gcf,[saveFolder '\GeoPosVelErr.eps'],...
        "ContentType","vector","Resolution",600)
end