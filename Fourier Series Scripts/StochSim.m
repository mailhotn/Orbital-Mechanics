%% Define Satellite & Noise Model
ic = [7000, 0.01, 55, 10, 10, 10];
primary = earth();
mu = primary.mu;
Sat = SingleSat(ic);
Prop = Propagator(Sat);

sigR = 5e-3;
sigV = 2e-5;
covEci = diag([sigR*ones(1,3),sigV*ones(1,3)]);
kMax = 5;
%% Prop & Contaminate
t = 0:100:86400;

[~, eciTrue] = Prop.PropEciJ2(t);
oeTrue = eci2oe(eciTrue);

noise = mvnrnd(zeros(6,1),covEci,length(t));
eciMeas = eciTrue + noise;
oeMeas = eci2oe(eciMeas);

%% Calculate Means
mOeTrue = osc2me(oeTrue);
mOeMeas = osc2me(oeMeas);
mOeFour = nan(size(oeTrue));

for iTime = 1:length(t)
    [~, lpeSpec] = LpeJ2Fourier(oeMeas(:,iTime),kMax);
    k = 1:kMax;
    M = oeMeas(6,iTime);
    a = oeMeas(1,iTime);
    nMo = sqrt(mu/a^3);
    trigVec = [cos(k*M)./k/nMo,sin(k*M)./k/nMo].';
    mOeFour(:,iTime) = oeMeas(:,iTime) - lpeSpec*trigVec;
end

%% Plot Stuff

figure(1)
plot(t,mOeMeas(1,:),t,mOeTrue(1,:),t,mOeFour(1,:))

figure(2)
plot(t,mOeMeas(2,:),t,mOeTrue(2,:))

figure(3)
plot(t,mOeMeas(3,:),t,mOeTrue(3,:))

figure(4)
plot(t,mOeMeas(4,:),t,mOeTrue(4,:))

figure(5)
plot(t,mOeMeas(5,:),t,mOeTrue(5,:))

figure(6)
plot(t,mOeMeas(6,:),t,mOeTrue(6,:))

