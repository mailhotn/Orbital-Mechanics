clear
%%
saveFolder = 'G:\My Drive\Doc Publications\2026 3 Body Fourier Paper\Figures';
% saveFolder = []; % comment to overwrite saved images

primary = Earth;
third = {Moon,Sun};
kMaxMol = 10;
kMaxGeo = 4;
dT = 1800; % 
nDay = 30;
t = 0:dT:86400*nDay;
%% Setup Molniya
T = 86164/2;
sma = ((T/2/pi)^2*primary.mu)^(1/3);
ic = [sma, 0.74, 63.4, 30, 270, 10].'; % Molniya
tOrb = 0:T:86164*nDay;

MolSat = SingleSat(ic,primary,third);
MolProp = Propagator(MolSat);
%% Propagate Molniya
% Prop using Good solution
tic
[~,X1] = MolProp.PropEci3B(t,'Stable');
oeN = eci2oe(X1.',[],primary,'me');
oeN(6,:) = unwrap(oeN(6,:)*pi/180)*180/pi - sqrt(primary.mu./oeN(1,:).^3).*t*180/pi;
stableTime = toc;
ra = oeN(1,:).*(1+oeN(2,:));
rp = oeN(1,:).*(1-oeN(2,:));
% oeN(1,:) = ra;
% oeN(2,:) = rp;
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
[~,oeF] = MolProp.PropOeFourier3B(t,kMaxMol);
oeF = oeF.';
fourTime = toc;
xF = oe2eci(oeF,Earth,'me');
% oeF(6,:) = me2ta(oeF(6,:),oeF(2,:));
oeF(6,:) = oeF(6,:) - sqrt(primary.mu./oeF(1,:).^3).*t*180/pi;
ra = oeF(1,:).*(1+oeF(2,:));
rp = oeF(1,:).*(1-oeF(2,:));
% oeF(1,:) = ra;
% oeF(2,:) = rp;
% lanF = oeF(5,:) + oeF(6,:)+ oeF(4,:);
% oeF(6,:) = lanF;

% Prop using singly averaged
tic
[~,oeA] = MolProp.PropOeMeanSingle3B(t);
oeA = oeA.';
singleTime = toc;
xA = oe2eci(oeA,Earth,'me');
oeA(6,:) = oeA(6,:) - sqrt(primary.mu./oeA(1,:).^3).*t*180/pi;
ra = oeA(1,:).*(1+oeA(2,:));
rp = oeA(1,:).*(1-oeA(2,:));
% oeA(1,:) = ra;
% oeA(2,:) = rp;
%% Import Stela Data
oeS = importdata('G:\My Drive\Doc Data\Stela 3B Sims\MolSimIcrf.mat').';
xS = oe2eci(oeS,Earth,'me');
xS(:,end) = nan(6,1);
oeS(6,:) = unwrap(oeS(6,:)*pi/180)*180/pi - sqrt(primary.mu./oeS(1,:).^3).*t*180/pi;
oeS(6,end) = nan; % weird outlier
ra = oeS(1,:).*(1+oeS(2,:));
rp = oeS(1,:).*(1-oeS(2,:));
% oeS(1,:) = ra;
% oeS(2,:) = rp;
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
dvF = [dF;vF];
dvA = [dA;vA];
dvS = [dS;vS];
rmsErrFMol = sqrt(trapz(t,errF.^2,2)/t(end))
rmsErrAMol = sqrt(trapz(t,errA.^2,2)/t(end))
rmsErrSMol = sqrt(trapz(t,errS.^2,2)/t(end))
molTimes = [stableTime,fourTime,singleTime]

%% Max error each orbit
maxErrF = nan(6,length(tOrb)-1);
maxErrA = nan(6,length(tOrb)-1);
maxErrS = nan(6,length(tOrb)-1);
maxDvF = nan(2,length(tOrb)-1);
maxDvA = nan(2,length(tOrb)-1);
maxDvS = nan(2,length(tOrb)-1);

for iOrb = 1:length(tOrb)-1
    orbInds = (t>=tOrb(iOrb) & t<tOrb(iOrb+1));
    maxErrF(:,iOrb) = max(errF(:,orbInds),[],2);
    maxErrA(:,iOrb) = max(errA(:,orbInds),[],2);
    maxErrS(:,iOrb) = max(errS(:,orbInds),[],2);
    maxDvF(:,iOrb) = max(dvF(:,orbInds),[],2);
    maxDvA(:,iOrb) = max(dvA(:,orbInds),[],2);
    maxDvS(:,iOrb) = max(dvS(:,orbInds),[],2);
end
%% Plot Molniya
% plot elements
lWidth = 1.5;
tD = t/86400;
iOrb = 1:length(tOrb)-1;
% figure("Units","centimeters","Position",[0,0,16,16])
% tiledlayout(3,2,"Padding","compact","Tilespacing","tight")
% nexttile
% plot(tD,oeF(1,:),tD,oeA(1,:),tD,oeS(1,:),':',tD,oeN(1,:),'--',LineWidth=lWidth)
% ylabel('$r_a \;[\rm{km}]$',Interpreter='latex',FontSize=12)
% xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
% nexttile
% plot(tD,oeF(2,:),tD,oeA(2,:),tD,oeS(2,:),':',tD,oeN(2,:),'--',LineWidth=lWidth)
% ylabel('$r_p \;[\rm{km}]$',Interpreter='latex',FontSize=12)
% xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
% legend('Fourier','Averaged','STELA','Numerical',Location='northwest')
% nexttile
% plot(tD,oeF(3,:),tD,oeA(3,:),tD,oeS(3,:),':',tD,oeN(3,:),'--',LineWidth=lWidth)
% ylabel('$i\; [\rm{deg}]$',Interpreter='latex',FontSize=12)
% xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
% nexttile
% plot(tD,oeF(4,:),tD,oeA(4,:),tD,oeS(4,:),':',tD,oeN(4,:),'--',LineWidth=lWidth)
% ylabel('$\Omega \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
% xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
% nexttile
% plot(tD,oeF(5,:),tD,oeA(5,:),tD,oeS(5,:),':',tD,oeN(5,:),'--',LineWidth=lWidth)
% ylabel('$\omega \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
% xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
% nexttile
% plot(tD,oeF(6,:),tD,oeA(6,:),tD,oeS(6,:),':',tD,oeN(6,:),'--',LineWidth=lWidth)
% ylabel('$M_0 \; [\rm{deg}]$',Interpreter='latex',FontSize=12)
% xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)

% if ~isempty(saveFolder)
% exportgraphics(gcf,[saveFolder '\MolniyaElements.eps'],...
%         "ContentType","vector","Resolution",600)
% end

% plot errors
figure("Units","centimeters","Position",[0,0,16,16])
tiledlayout(3,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(iOrb,maxErrF(1,:),'o',iOrb,maxErrA(1,:),'x',iOrb,maxErrS(1,:),'d',LineWidth=lWidth)
ylabel('$\delta a \;[\rm{km}]$',Interpreter='latex',FontSize=12)
xlabel('$\rm{Orbit}$',Interpreter='latex',FontSize=12)
nexttile
plot(iOrb,maxErrF(2,:),'o',iOrb,maxErrA(2,:),'x',iOrb,maxErrS(2,:),'d',LineWidth=lWidth)
ylabel('$\delta e$',Interpreter='latex',FontSize=12)
xlabel('$\rm{Orbit}$',Interpreter='latex',FontSize=12)
legend('Fourier','Averaged','STELA',Location='northwest')
nexttile
plot(iOrb,maxErrF(3,:),'o',iOrb,maxErrA(3,:),'x',iOrb,maxErrS(3,:),'d',LineWidth=lWidth)
ylabel('$\delta i \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$\rm{Orbit}$',Interpreter='latex',FontSize=12)
nexttile
plot(iOrb,maxErrF(4,:),'o',iOrb,maxErrA(4,:),'x',iOrb,maxErrS(4,:),'d',LineWidth=lWidth)
ylabel('$\delta \Omega \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$\rm{Orbit}$',Interpreter='latex',FontSize=12)
nexttile
plot(iOrb,maxErrF(5,:),'o',iOrb,maxErrA(5,:),'x',iOrb,maxErrS(5,:),'d',LineWidth=lWidth)
ylabel('$\delta \omega \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$\rm{Orbit}$',Interpreter='latex',FontSize=12)
nexttile
plot(iOrb,maxErrF(6,:),'o',iOrb,maxErrA(6,:),'x',iOrb,maxErrS(6,:),'d',LineWidth=lWidth)
ylabel('$\delta M_0 \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$\rm{Orbit}$',Interpreter='latex',FontSize=12)

if ~isempty(saveFolder)
exportgraphics(gcf,[saveFolder '\MolniyaErrors.eps'],...
        "ContentType","vector","Resolution",600)
end

% plot dist vel errors
figure("Units","centimeters","Position",[0,0,16,6])
tiledlayout(1,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(iOrb,maxDvF(1,:),'o',iOrb,maxDvA(1,:),'x',iOrb,maxDvS(1,:),'d',LineWidth=lWidth)
ylabel('$\delta d \;[\rm{km}]$',Interpreter='latex',FontSize=12)
xlabel('$\rm{Orbit}$',Interpreter='latex',FontSize=12)
nexttile
plot(iOrb,maxDvF(2,:),'o',iOrb,maxDvA(2,:),'x',iOrb,maxDvS(2,:),'d',LineWidth=lWidth)
ylabel('$\delta v \;\left[\rm{km}\,\rm{s}^{-1}\right]$',Interpreter='latex',FontSize=12)
xlabel('$\rm{Orbit}$',Interpreter='latex',FontSize=12)
legend('Fourier','Averaged','STELA',Location='northwest')

if ~isempty(saveFolder)
exportgraphics(gcf,[saveFolder '\MolniyaPosVelErr.eps'],...
        "ContentType","vector","Resolution",600)
end
%% Setup GEO
T = 86164;
sma = ((T/2/pi)^2*primary.mu)^(1/3);
ic = [sma, 0.001, 5, 30, 30, 10].'; % GEO
tOrb = 0:T:86400*nDay;

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
ra = oeN(1,:).*(1+oeN(2,:));
rp = oeN(1,:).*(1-oeN(2,:));
% oeN(1,:) = ra;
% oeN(2,:) = rp;

% Prop using Fourier solution
tic
[~,oeF] = GeoProp.PropOeFourier3B(t,kMaxGeo);
oeF = oeF.';
fourTime = toc;
xF = oe2eci(oeF,Earth,'me');
oeF(6,:) = oeF(6,:) - sqrt(primary.mu./oeF(1,:).^3).*t*180/pi;
lanF = oeF(5,:) + oeF(6,:)+ oeF(4,:);
oeF(4,:) = lanF;
ra = oeF(1,:).*(1+oeF(2,:));
rp = oeF(1,:).*(1-oeF(2,:));
% oeF(1,:) = ra;
% oeF(2,:) = rp;

% Prop using singly averaged
tic
[~,oeA] = GeoProp.PropOeMeanSingle3B(t);
oeA = oeA.';
singleTime = toc;
xA = oe2eci(oeA,Earth,'me');
oeA(6,:) = oeA(6,:) - sqrt(primary.mu./oeA(1,:).^3).*t*180/pi;
lanA =  oeA(5,:) + oeA(6,:) + oeA(4,:);
oeA(4,:) = lanA;
ra = oeA(1,:).*(1+oeA(2,:));
rp = oeA(1,:).*(1-oeA(2,:));
% oeA(1,:) = ra;
% oeA(2,:) = rp;
%% Import Stela Data
oeS = importdata('G:\My Drive\Doc Data\Stela 3B Sims\GeoSimIcrf.mat').';
xS = oe2eci(oeS,Earth,'me');
xS(:,end) = nan(6,1);
oeS(6,:) = unwrap(oeS(6,:)*pi/180)*180/pi - sqrt(primary.mu./oeS(1,:).^3).*t*180/pi;
lanS = oeS(5,:) + oeS(6,:)+oeS(4,:);
lanS(end) = nan; % weird outlier
oeS(4,:) = lanS;
ra = oeS(1,:).*(1+oeS(2,:));
rp = oeS(1,:).*(1-oeS(2,:));
% oeS(1,:) = ra;
% oeS(2,:) = rp;
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
dvF = [dF;vF];
dvA = [dA;vA];
dvS = [dS;vS];
rmsErrFGeo = sqrt(trapz(t,errF.^2,2)/t(end))
rmsErrAGeo = sqrt(trapz(t,errA.^2,2)/t(end))
rmsErrSGeo = sqrt(trapz(t,errS.^2,2)/t(end))
geoTimes = [stableTime,fourTime,singleTime]
%% Max error each orbit
maxErrF = nan(4,length(tOrb)-1);
maxErrA = nan(4,length(tOrb)-1);
maxErrS = nan(4,length(tOrb)-1);
maxDvF = nan(2,length(tOrb)-1);
maxDvA = nan(2,length(tOrb)-1);
maxDvS = nan(2,length(tOrb)-1);

for iOrb = 1:length(tOrb)-1
    orbInds = (t>=tOrb(iOrb) & t<tOrb(iOrb+1));
    maxErrF(:,iOrb) = max(errF(1:4,orbInds),[],2);
    maxErrA(:,iOrb) = max(errA(1:4,orbInds),[],2);
    maxErrS(:,iOrb) = max(errS(1:4,orbInds),[],2);
    maxDvF(:,iOrb) = max(dvF(:,orbInds),[],2);
    maxDvA(:,iOrb) = max(dvA(:,orbInds),[],2);
    maxDvS(:,iOrb) = max(dvS(:,orbInds),[],2);
end
%% Plot GEO
% plot elements
lWidth = 1.5;
tD = t/86400;
iOrb = 1:length(tOrb)-1;
% figure("Units","centimeters","Position",[0,0,16,12])
% tiledlayout(2,2,"Padding","compact","Tilespacing","tight")
% nexttile
% plot(tD,oeF(1,:),tD,oeA(1,:),tD,oeS(1,:),':',tD,oeN(1,:),'--',LineWidth=lWidth)
% ylabel('$r_a \;[\rm{km}]$',Interpreter='latex',FontSize=12)
% xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
% nexttile
% plot(tD,oeF(2,:),tD,oeA(2,:),tD,oeS(2,:),':',tD,oeN(2,:),'--',LineWidth=lWidth)
% ylabel('$r_p \;[\rm{km}]$',Interpreter='latex',FontSize=12)
% xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
% nexttile
% plot(tD,oeF(3,:),tD,oeA(3,:),tD,oeS(3,:),':',tD,oeN(3,:),'--',LineWidth=lWidth)
% ylabel('$i\; [\rm{deg}]$',Interpreter='latex',FontSize=12)
% xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
% nexttile
% plot(tD,lanF,tD,lanA,tD,lanS,':',tD,lanN,'--',LineWidth=lWidth)
% ylabel('$L \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
% xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
% legend('Fourier','Averaged','STELA','Numerical',Location='southwest')
% 
% if ~isempty(saveFolder)
% exportgraphics(gcf,[saveFolder '\GeoElements.eps'],...
%         "ContentType","vector","Resolution",600)
% end

% plot errors
figure("Units","centimeters","Position",[0,0,16,12])
tiledlayout(2,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(iOrb,maxErrF(1,:),'o',iOrb,maxErrA(1,:),'x',iOrb,maxErrS(1,:),'d',LineWidth=lWidth)
ylabel('$\delta a \;[\rm{km}]$',Interpreter='latex',FontSize=12)
xlabel('$\rm{Orbit}$',Interpreter='latex',FontSize=12)
nexttile
plot(iOrb,maxErrF(2,:),'o',iOrb,maxErrA(2,:),'x',iOrb,maxErrS(2,:),'d',LineWidth=lWidth)
ylabel('$\delta e $',Interpreter='latex',FontSize=12)
xlabel('$\rm{Orbit}$',Interpreter='latex',FontSize=12)
legend('Fourier','Averaged','STELA',Location='northwest')
nexttile
plot(iOrb,maxErrF(3,:),'o',iOrb,maxErrA(3,:),'x',iOrb,maxErrS(3,:),'d',LineWidth=lWidth)
ylabel('$\delta i \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$\rm{Orbit}$',Interpreter='latex',FontSize=12)
nexttile
plot(iOrb,maxErrF(4,:),'o',iOrb,maxErrA(4,:),'x',iOrb,maxErrS(4,:),'d',LineWidth=lWidth)
ylabel('$\delta L \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$\rm{Orbit}$',Interpreter='latex',FontSize=12)

if ~isempty(saveFolder)
exportgraphics(gcf,[saveFolder '\GeoErrors.eps'],...
        "ContentType","vector","Resolution",600)
end

% plot dist vel errors
figure("Units","centimeters","Position",[0,0,16,6])
tiledlayout(1,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(iOrb,maxDvF(1,:),'o',iOrb,maxDvA(1,:),'x',iOrb,maxDvS(1,:),'d',LineWidth=lWidth)
ylabel('$\delta d \;[\rm{km}]$',Interpreter='latex',FontSize=12)
xlabel('$\rm{Orbit}$',Interpreter='latex',FontSize=12)
nexttile
plot(iOrb,maxDvF(2,:),'o',iOrb,maxDvA(2,:),'x',iOrb,maxDvS(2,:),'d',LineWidth=lWidth)
ylabel('$\delta v \;\left[\rm{km}\,\rm{s}^{-1}\right]$',Interpreter='latex',FontSize=12)
xlabel('$\rm{Orbit}$',Interpreter='latex',FontSize=12)
legend('Fourier','Averaged','STELA',Location='northwest')

if ~isempty(saveFolder)
exportgraphics(gcf,[saveFolder '\GeoPosVelErr.eps'],...
        "ContentType","vector","Resolution",600)
end