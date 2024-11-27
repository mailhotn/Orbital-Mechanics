clear
oeLead = [7000, 0.005, 5, 5, 5, 5].'; % init oe
dOe = [0,0,0,0,0,1].'; % mean element difference for follower
oeMat = [oeLead,oeLead + dOe];
kMax = 5;

Form = GeneralFormation(oeLead,dOe);
Prop = Propagator(Form);
% leadSat = SingleSat(oeLead);
% leadProp = Propagator(leadSat);
% followSat = SingleSat(oeLead+dOe);
% followProp = Propagator(followSat);

primary = earth();
mu = primary.mu;

nDay = 1;
t = 0:100:86400*nDay;

%% Sim Difference Prop
% Prop Formation
[~,oeC] = Prop.PropOeOsc3(t);
oeC = oeC.';
[~,oeM] = Prop.PropOeMeanShort(t);
oeB = reshape(me2oscSP(reshape(oeM.',6,length(t)*2)),12,length(t));
[~,oeF] = Prop.PropOeFourier2Ord(t,kMax);
oeF = oeF.';

% Get relative ECI vector
dEciC = oe2eci(oeC(7:12,:))-oe2eci(oeC(1:6,:));
dEciB = oe2eci(oeB(7:12,:))-oe2eci(oeB(1:6,:));
dEciF = oe2eci(oeF(7:12,:))-oe2eci(oeF(1:6,:));
% Distance
distC = vecnorm(dEciC(1:3,:));
distB = vecnorm(dEciB(1:3,:));
distF = vecnorm(dEciF(1:3,:));
% Velocity
relVC = vecnorm(dEciC(4:6,:));
relVB = vecnorm(dEciB(4:6,:));
relVF = vecnorm(dEciF(4:6,:));
% Differential elements
dOeC = oeC(7:12,:)-oeC(1:6,:);
dOeB = oeB(7:12,:)-oeB(1:6,:);
dOeF = oeF(7:12,:)-oeF(1:6,:);
% Diff OE errors, normalized by nominal a,e, radians
errOeB = dOeB-dOeC;
errOeB = [errOeB(1:2,:)./oeLead(1:2); errOeB(3:6,:)*pi/180];
errOeF = dOeF-dOeC;
errOeF = [errOeF(1:2,:)./oeLead(1:2); errOeF(3:6,:)*pi/180];
% Distance & velocity errors
errDB = distB-distC;
errDF = distF-distC;

errVB = relVB-relVC;
errVF = relVF-relVC;

trapz(t,abs(errOeB),2)/t(end)
trapz(t,abs(errOeF),2)/t(end)

trapz(t,abs(errDB))/t(end)*1000
trapz(t,abs(errDF))/t(end)*1000

trapz(t,abs(errVB))/t(end)*1000*100
trapz(t,abs(errVF))/t(end)*1000*100
%% Plot
% figure(1)
% plot(t/86400,distF)
% 
% figure(2)
% plot(t/86400,distB)

% Relative Distance & Velocity
figure(1)
plot(t/3600,errDB,t/3600,errDF)
legend('Brouwer','Fourier')
xlabel('$Time \left[\rm{hr}\right]$',Interpreter='latex',FontSize=12)
xlim([0,nDay*24])
xticks(linspace(0,nDay*24,5*nDay))
ylabel('$\delta r \ Error \left[\rm{km}\right]$',Interpreter='latex',FontSize=12)

figure(2)
plot(t/3600,errVB,t/3600,errVF)
legend('Brouwer','Fourier')
xlabel('$Time \left[\rm{hr}\right]$',Interpreter='latex',FontSize=12)
xlim([0,nDay*24])
xticks(linspace(0,nDay*24,5*nDay))
ylabel('$\delta v \ Error \left[\rm{\frac{km}{s}}\right]$',Interpreter='latex',FontSize=12)

figure(3)
plot(t/3600,dEciB(1,:)-dEciC(1,:),t/3600,dEciF(1,:)-dEciC(1,:))
legend('Brouwer','Fourier')
xlabel('$Time \left[\rm{hr}\right]$',Interpreter='latex',FontSize=12)
xlim([0,nDay*24])
xticks(linspace(0,nDay*24,5*nDay))
ylabel('$\delta x \ Error \left[\rm{km}\right]$',Interpreter='latex',FontSize=12)
figure(4)
plot(t/3600,dEciB(2,:)-dEciC(2,:),t/3600,dEciF(2,:)-dEciC(2,:))
legend('Brouwer','Fourier')
xlabel('$Time \left[\rm{hr}\right]$',Interpreter='latex',FontSize=12)
xlim([0,nDay*24])
xticks(linspace(0,nDay*24,5*nDay))
ylabel('$\delta y \ Error \left[\rm{km}\right]$',Interpreter='latex',FontSize=12)
figure(5)
plot(t/3600,dEciB(3,:)-dEciC(3,:),t/3600,dEciF(3,:)-dEciC(3,:))
legend('Brouwer','Fourier')
xlabel('$Time \left[\rm{hr}\right]$',Interpreter='latex',FontSize=12)
xlim([0,nDay*24])
xticks(linspace(0,nDay*24,5*nDay))
ylabel('$\delta z \ Error \left[\rm{km}\right]$',Interpreter='latex',FontSize=12)

% orbital Elements
figure(11)
plot(t/3600,errOeB(1,:),t/3600,errOeF(1,:))
legend('Brouwer','Fourier')
xlabel('$Time \left[\rm{hr}\right]$',Interpreter='latex',FontSize=12)
xlim([0,nDay*24])
xticks(linspace(0,nDay*24,5*nDay))
ylabel('$ \delta a \ Error$',Interpreter='latex',FontSize=12)

figure(12)
plot(t/3600,errOeB(2,:),t/3600,errOeF(2,:))
legend('Brouwer','Fourier')
xlabel('$Time \left[\rm{hr}\right]$',Interpreter='latex',FontSize=12)
xlim([0,nDay*24])
xticks(linspace(0,nDay*24,5*nDay))
ylabel('$ \delta e \ Error$',Interpreter='latex',FontSize=12)

figure(13)
plot(t/3600,errOeB(3,:),t/3600,errOeF(3,:))
legend('Brouwer','Fourier')
xlabel('$Time \left[\rm{hr}\right]$',Interpreter='latex',FontSize=12)
xlim([0,nDay*24])
xticks(linspace(0,nDay*24,5*nDay))
ylabel('$ \delta i \ Error \left[\rm{rad}\right]$',Interpreter='latex',FontSize=12)

figure(14)
plot(t/3600,errOeB(4,:),t/3600,errOeF(4,:))
legend('Brouwer','Fourier')
xlabel('$Time \left[\rm{hr}\right]$',Interpreter='latex',FontSize=12)
xlim([0,nDay*24])
xticks(linspace(0,nDay*24,5*nDay))
ylabel('$ \delta \Omega \ Error \left[\rm{rad}\right]$',Interpreter='latex',FontSize=12)

figure(15)
plot(t/3600,errOeB(5,:),t/3600,errOeF(5,:))
legend('Brouwer','Fourier')
xlabel('$Time \left[\rm{hr}\right]$',Interpreter='latex',FontSize=12)
xlim([0,nDay*24])
xticks(linspace(0,nDay*24,5*nDay))
ylabel('$ \delta \omega \ Error \left[\rm{rad}\right]$',Interpreter='latex',FontSize=12)

figure(16)
plot(t/3600,errOeB(6,:),t/3600,errOeF(6,:))
legend('Brouwer','Fourier')
xlabel('$Time \left[\rm{hr}\right]$',Interpreter='latex',FontSize=12)
xlim([0,nDay*24])
xticks(linspace(0,nDay*24,5*nDay))
ylabel('$ \delta M \ Error \left[\rm{rad}\right]$',Interpreter='latex',FontSize=12)

%% Init Difference - Doesn't work
% % Calculate osc elements w/ Brouwer & w/ Fourier
% oeMatB = me2osc(oeMat);
% 
% % Fourier
% [~,lpeSpec] = LpeJ2Fourier(oeMat,kMax);
% M = oeMat(6,:);
% nomSma = oeLead(1);
% nMo = sqrt(mu./nomSma.^3);
% k = (1:kMax).';
% Ck = -cos(k*M)./k./nMo;
% Sk = sin(k*M)./k./nMo;
% 
% % Fix M
% manVar =  sum(Sk.*squeeze(lpeSpec(11,:,:)).',1) + ...
%     sum(Ck.*squeeze(lpeSpec(12,:,:)).',1);
% M2 = oeMat(6,:) + manVar; % Fixed osc M
% Ck = -cos(k*M2)./k./nMo;
% Sk = sin(k*M2)./k./nMo;
% % Calculate other Elements
% smaVar =  sum(Sk.*squeeze(lpeSpec(1,:,:)).',1) + ...
%     sum(Ck.*squeeze(lpeSpec(2,:,:)).',1);
% eccVar =  sum(Sk.*squeeze(lpeSpec(3,:,:)).',1) + ...
%     sum(Ck.*squeeze(lpeSpec(4,:,:)).',1);
% incVar =  sum(Sk.*squeeze(lpeSpec(5,:,:)).',1) + ...
%     sum(Ck.*squeeze(lpeSpec(6,:,:)).',1);
% ranVar =  sum(Sk.*squeeze(lpeSpec(7,:,:)).',1) + ...
%     sum(Ck.*squeeze(lpeSpec(8,:,:)).',1);
% aopVar =  sum(Sk.*squeeze(lpeSpec(9,:,:)).',1) + ...
%     sum(Ck.*squeeze(lpeSpec(10,:,:)).',1);
% manVar =  sum(Sk.*squeeze(lpeSpec(11,:,:)).',1) + ...
%     sum(Ck.*squeeze(lpeSpec(12,:,:)).',1);
% 
% oeMatF = oeMat + [smaVar;eccVar;incVar;ranVar;aopVar;manVar];
% 
% %% Prop both scenarios with OeOsc3
% bForm = GeneralFormation(oeMatB(:,1),oeMatB(:,2)-oeMatB(:,1),false);
% bProp = Propagator(bForm);
% [~,oeB] = bProp.PropOeOsc3(t);
% oeB = me2ta(reshape(oeB.',6,length(t)*2));
% xB = reshape(oe2eci(oeB),12,length(t));
% distB = vecnorm(xB(7:9,:) - xB(1:3,:));
% 
% fForm = GeneralFormation(oeMatF(:,1),oeMatF(:,2)-oeMatF(:,1),false);
% fProp = Propagator(fForm);
% [~,oeF] = fProp.PropOeOsc3(t);
% oeF = me2ta(reshape(oeF.',6,length(t)*2));
% xF = reshape(oe2eci(oeF),12,length(t));
% distF = vecnorm(xF(7:9,:) - xF(1:3,:)); 

