clear
oeLead = [7000, 0.01, 5, 0, 0, 0].'; % init oe
dOe = [0,0,0,0,0,1].'; % mean element difference for follower
oeMat = [oeLead,oeLead + dOe];
kMax = 5;

leadSat = SingleSat(oeLead);
leadProp = Propagator(leadSat);

primary = earth();
mu = primary.mu;

t = 0:100:86400*5;

%% Sim Difference Prop
[~,oeC] = leadProp.PropOeOsc3(t);
oeC = oeC.';
[~,oeM] = leadProp.PropOeMeanShort(t);
oeB = me2oscSP(oeM.');
[~,oeF] = leadProp.PropOeFourier2Ord(t,5);
oeF = oeF.';
errB = abs(oeC-oeB);
errF = abs(oeC-oeF);

%% Init Difference
% Calculate osc elements w/ Brouwer & w/ Fourier
oeMatB = me2osc(oeMat);

% Fourier
[~,lpeSpec] = LpeJ2Fourier(oeMat,kMax);
M = oeMat(6,:);
nomSma = oeLead(1);
nMo = sqrt(mu./nomSma.^3);
k = (1:kMax).';
Ck = -cos(k*M)./k./nMo;
Sk = sin(k*M)./k./nMo;

% Fix M
manVar =  sum(Sk.*squeeze(lpeSpec(11,:,:)).',1) + ...
    sum(Ck.*squeeze(lpeSpec(12,:,:)).',1);
M2 = oeMat(6,:) + manVar; % Fixed osc M
Ck = -cos(k*M2)./k./nMo;
Sk = sin(k*M2)./k./nMo;
% Calculate other Elements
smaVar =  sum(Sk.*squeeze(lpeSpec(1,:,:)).',1) + ...
    sum(Ck.*squeeze(lpeSpec(2,:,:)).',1);
eccVar =  sum(Sk.*squeeze(lpeSpec(3,:,:)).',1) + ...
    sum(Ck.*squeeze(lpeSpec(4,:,:)).',1);
incVar =  sum(Sk.*squeeze(lpeSpec(5,:,:)).',1) + ...
    sum(Ck.*squeeze(lpeSpec(6,:,:)).',1);
ranVar =  sum(Sk.*squeeze(lpeSpec(7,:,:)).',1) + ...
    sum(Ck.*squeeze(lpeSpec(8,:,:)).',1);
aopVar =  sum(Sk.*squeeze(lpeSpec(9,:,:)).',1) + ...
    sum(Ck.*squeeze(lpeSpec(10,:,:)).',1);
manVar =  sum(Sk.*squeeze(lpeSpec(11,:,:)).',1) + ...
    sum(Ck.*squeeze(lpeSpec(12,:,:)).',1);

oeMatF = oeMat + [smaVar;eccVar;incVar;ranVar;aopVar;manVar];

%% Prop both scenarios with OeOsc3
bForm = GeneralFormation(oeMatB(:,1),oeMatB(:,2)-oeMatB(:,1),false);
bProp = Propagator(bForm);
[~,oeB] = bProp.PropOeOsc3(t);
oeB = me2ta(reshape(oeB.',6,length(t)*2));
xB = reshape(oe2eci(oeB),12,length(t));
distB = vecnorm(xB(7:9,:) - xB(1:3,:));

fForm = GeneralFormation(oeMatF(:,1),oeMatF(:,2)-oeMatF(:,1),false);
fProp = Propagator(fForm);
[~,oeF] = fProp.PropOeOsc3(t);
oeF = me2ta(reshape(oeF.',6,length(t)*2));
xF = reshape(oe2eci(oeF),12,length(t));
distF = vecnorm(xF(7:9,:) - xF(1:3,:)); 

%% Plot
figure(1)
plot(t/86400,distF)

figure(2)
plot(t/86400,distB)