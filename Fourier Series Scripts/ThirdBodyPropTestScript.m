clear
%%
saveFolder = 'G:\My Drive\Doc Publications\2026 3 Body Fourier Paper\Figures';
saveFolder = [];
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
% oeF(6,:) = me2ta(oeF(6,:),oeF(2,:));
oeF(6,:) = oeF(6,:) - sqrt(primary.mu./oeF(1,:).^3).*t*180/pi;
% lanF = oeF(5,:) + oeF(6,:)+ oeF(4,:);
% oeF(6,:) = lanF;

% Prop using singly averaged
tic
[~,oeS] = MolProp.PropOeMeanSingle3B(t);
oeS = oeS.';
singleTime = toc;
oeS(6,:) = oeS(6,:) - sqrt(primary.mu./oeS(1,:).^3).*t*180/pi;
% lanS =  oeS(5,:) + oeS(6,:) + oeS(4,:);
% oeS(6,:) = lan2;
%% Calculate relative errors - Molniya
errS = abs(oeS-oeN)./abs(oeN);
% err3 = abs(oeN)*0;
errF = abs(oeF-oeN)./abs(oeN);
meanErrFMol = trapz(t,errF,2)/t(end)
meanErrSMol = trapz(t,errS,2)/t(end)
molTimes = [stableTime,fourTime,singleTime]
%% Plot Molniya
% plot elements
lWidth = 1.5;
tD = t/86400;
figure("Units","centimeters","Position",[0,0,16,16])
tiledlayout(3,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(tD,oeF(1,:),tD,oeS(1,:),tD,oeN(1,:),'--',LineWidth=lWidth)
ylabel('$a \;[\rm{km}]$',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,oeF(2,:),tD,oeS(2,:),tD,oeN(2,:),'--',LineWidth=lWidth)
ylabel('$e$',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
legend('Fourier','Singly Averaged','Numerical',Location='southeast')
nexttile
plot(tD,oeF(3,:),tD,oeS(3,:),tD,oeN(3,:),'--',LineWidth=lWidth)
ylabel('$i\; [\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,oeF(4,:),tD,oeS(4,:),tD,oeN(4,:),'--',LineWidth=lWidth)
ylabel('$\Omega \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,oeF(5,:),tD,oeS(5,:),tD,oeN(5,:),'--',LineWidth=lWidth)
ylabel('$\omega \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,oeF(6,:),tD,oeS(6,:),tD,oeN(6,:),'--',LineWidth=lWidth)
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
plot(tD,errF(1,:),tD,errS(1,:),'--',LineWidth=lWidth)
ylabel('$\delta a \;[\rm{km}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,errF(2,:),tD,errS(2,:),'--',LineWidth=lWidth)
ylabel('$\delta e $',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
legend('Fourier','Singly Averaged')
nexttile
plot(tD,errF(3,:),tD,errS(3,:),'--',LineWidth=lWidth)
ylabel('$\delta i \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,errF(4,:),tD,errS(4,:),'--',LineWidth=lWidth)
ylabel('$\delta \Omega \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,errF(5,:),tD,errS(5,:),'--',LineWidth=lWidth)
ylabel('$\delta \omega \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,errF(6,:),tD,errS(6,:),'--',LineWidth=lWidth)
ylabel('$\delta M_0 \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)

if ~isempty(saveFolder)
exportgraphics(gcf,[saveFolder '\MolniyaErrors.eps'],...
        "ContentType","vector","Resolution",600)
end
% figure(3)
% bar([stableTime,fourTime,singleTime])

% figure(4)
% plot(tD,lan1,tD,lan2,tD,lanF,'--',LineWidth=lWidth)
% ylabel('L - nt [\rm{deg}]')
% xlabel('t [\rm{day}]')
% legend('Numerical - Stable','Numerical - Taylor','Fourier')
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
oeF(6,:) = oeF(6,:) - sqrt(primary.mu./oeF(1,:).^3).*t*180/pi;
lanF = oeF(5,:) + oeF(6,:)+ oeF(4,:);
oeF(4,:) = lanF;

% Prop using singly averaged
tic
[~,oeS] = GeoProp.PropOeMeanSingle3B(t);
oeS = oeS.';
singleTime = toc;
oeS(6,:) = oeS(6,:) - sqrt(primary.mu./oeS(1,:).^3).*t*180/pi;
lanS =  oeS(5,:) + oeS(6,:) + oeS(4,:);
oeS(4,:) = lanS;
%% Calculate relative errors - GEO
errS = abs(oeS-oeN)./abs(oeN);
err3 = abs(oeN)*0;
errF = abs(oeF-oeN)./abs(oeN);
meanErrFGeo = trapz(t,errF(1:4,:),2)/t(end)
meanErrSGeo = trapz(t,errS(1:4,:),2)/t(end)
geoTimes = [stableTime,fourTime,singleTime]
%% Plot GEO
% plot elements
lWidth = 1.5;
tD = t/86400;
figure("Units","centimeters","Position",[0,0,16,12])
tiledlayout(2,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(tD,oeF(1,:),tD,oeS(1,:),tD,oeN(1,:),'--',LineWidth=lWidth)
ylabel('$a \;[\rm{km}]$',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,oeF(2,:),tD,oeS(2,:),tD,oeN(2,:),'--',LineWidth=lWidth)
ylabel('$e$',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
legend('Fourier','Singly Averaged','Numerical',Location='southeast')
nexttile
plot(tD,oeF(3,:),tD,oeS(3,:),tD,oeN(3,:),'--',LineWidth=lWidth)
ylabel('$i\; [\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,lanF,tD,lanS,tD,lanN,'--',LineWidth=lWidth)
ylabel('$L \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)

if ~isempty(saveFolder)
exportgraphics(gcf,[saveFolder '\GeoElements.eps'],...
        "ContentType","vector","Resolution",600)
end

% plot errors
figure("Units","centimeters","Position",[0,0,16,12])
tiledlayout(2,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(tD,errF(1,:),tD,errS(1,:),'--',LineWidth=lWidth)
ylabel('$\delta a \;[\rm{km}]$',Interpreter='latex',FontSize=12)
xlabel('$t\; [\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,errF(2,:),tD,errS(2,:),'--',LineWidth=lWidth)
ylabel('$\delta e $',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
legend('Fourier','Singly Averaged')
nexttile
plot(tD,errF(3,:),tD,errS(3,:),'--',LineWidth=lWidth)
ylabel('$\delta i \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)
nexttile
plot(tD,errF(4,:),tD,errS(4,:),'--',LineWidth=lWidth)
ylabel('$\delta L \;[\rm{deg}]$',Interpreter='latex',FontSize=12)
xlabel('$t \;[\rm{day}]$',Interpreter='latex',FontSize=12)

if ~isempty(saveFolder)
exportgraphics(gcf,[saveFolder '\GeoErrors.eps'],...
        "ContentType","vector","Resolution",600)
end