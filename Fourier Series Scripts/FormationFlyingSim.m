clear
oeLead = [7000, 0.01, 5, 10, 10, 5].'; % init oe
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

t = 0:100:86400*5;

%% Sim Difference Prop
% Prop Formation
[~,oeC] = Prop.PropOeOsc3(t);
oeC = oeC.';
[~,oeM] = Prop.PropOeMeanShort(t);
oeB = reshape(me2oscSP(reshape(oeM.',6,length(t)*2)),12,length(t));
[~,oeF] = Prop.PropOeFourier2Ord(t,kMax);
oeF = oeF.';

% Get relative position
dEciC = oe2eci(oeC(7:12,:))-oe2eci(oeC(1:6,:));
dEciB = oe2eci(oeB(7:12,:))-oe2eci(oeB(1:6,:));
dEciF = oe2eci(oeF(7:12,:))-oe2eci(oeF(1:6,:));

distC = vecnorm(dEciC(1:3,:));
distB = vecnorm(dEciB(1:3,:));
distF = vecnorm(dEciF(1:3,:));

dOeC = oeC(7:12,:)-oeC(1:6,:);
dOeB = oeB(7:12,:)-oeB(1:6,:);
dOeF = oeF(7:12,:)-oeF(1:6,:);

errOeB = dOeB-dOeC;
errOeF = dOeF-dOeC;

errB = distB-distC;
errF = distF-distC;

%% Plot
figure(1)
plot(t/86400,distF)

figure(2)
plot(t/86400,distB)

figure(3)
plot(t/86400,errB,t/86400,errF)

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

