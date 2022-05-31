clear
nSats = 1000; % different orbits
kMaxMax = 7;
nOrb = 5; % # orbits for test

aErrMat = nan(nSats,kMaxMax+1);
eErrMat = nan(nSats,kMaxMax+1);
iErrMat = nan(nSats,kMaxMax+1);
raanErrMat = nan(nSats,kMaxMax+1);
aopErrMat = nan(nSats,kMaxMax+1);
MErrMat = nan(nSats,kMaxMax+1);

for iSat = 1:nSats
    % random orbit
    r = rand(6,1);
    oe = [7500+r(1)*3000;
        0.01+0.1*r(2);
        63.4+0.5*(2*r(3)-1);
        r(4)*360;
        r(5)*360;
        r(6)*360].';
    
    % Define Sat & Prop
    Sat = SingleSat(oe,earth());
    Prop = Propagator(Sat);
    T = 2*pi*sqrt(oe(1)^3/Sat.primary.mu);
    t = 0:100:nOrb*T;
    
    % Prop Numerical
    [~,oeC] = Prop.PropOeOsc(t);
    oeC = oeC.';
    
    % Prop Brouwer & Get Error
    [~,OeM] = Prop.PropOeMeanFast(t);
    oeB = me2osc(OeM.');
    
    errB = abs(oeC-oeB);
    errB = [errB(1,:)/oe(1);errB(2,:)/oe(2);errB(3:end,:)*pi/180];
    intErrB = trapz(t.',errB,2)/t(end);
    aErrMat(iSat,1) = intErrB(1);
    eErrMat(iSat,1) = intErrB(2);
    iErrMat(iSat,1) = intErrB(3);
    raanErrMat(iSat,1) = intErrB(4);
    aopErrMat(iSat,1) = intErrB(5);
    MErrMat(iSat,1) = intErrB(6);
    
    for kMax = 1:kMaxMax
        % Prop Fourier & get error
        [~,oeF] = Prop.PropOeFourier2(t,kMax);
        errF = abs(oeC-oeF);
        errF = [errF(1,:)/oe(1);errF(2,:)/oe(2);errF(3:end,:)*pi/180];
        intErrF = trapz(t.',errF,2)/t(end);
        aErrMat(iSat,kMax+1) = intErrF(1);
        eErrMat(iSat,kMax+1) = intErrF(2);
        iErrMat(iSat,kMax+1) = intErrF(3);
        raanErrMat(iSat,kMax+1) = intErrF(4);
        aopErrMat(iSat,kMax+1) = intErrF(5);
        MErrMat(iSat,kMax+1) = intErrF(6);
    end
end

%% Plot - Errorbar stuff
k = 1:kMaxMax;
% aStd = std(aErrMat);
% eStd = std(eErrMat);
% iStd = std(iErrMat);
% raanStd = std(raanErrMat);
% aopStd = std(aopErrMat);
% MStd = std(MErrMat);



aErrM = mean(aErrMat);
eErrM = mean(eErrMat);
iErrM = mean(iErrMat);
raanErrM = mean(raanErrMat);
aopErrM = mean(aopErrMat);
MErrM = mean(MErrMat);

figure(1)
clf
errorbar(k,aErrM(2:end),aErrM(2:end)-prctile(aErrMat(2:end),20),prctile(aErrMat(2:end),80)-aErrM(2:end),'o','linewidth',2)
hold on
plot(0:(kMaxMax+1),aErrM(1)*ones(1,length(k)+2),'--','linewidth',2)
xticks(1:kMaxMax)
xlim([0,kMaxMax+1])
grid on
ylabel('$\textrm{Average error}$','interpreter','latex','fontsize',18)
xlabel('$k_{Max}$','interpreter','latex','fontsize',18)
legend('Fourier','Brouwer','fontsize',12)
hold off

figure(2)
clf
errorbar(k,eErrM(2:end),eErrM(2:end)-prctile(eErrMat(2:end),20),prctile(eErrMat(2:end),80)-eErrM(2:end),'o','linewidth',2)
hold on
plot(0:(kMaxMax+1),eErrM(1)*ones(1,length(k)+2),'--','linewidth',2)
xticks(1:kMaxMax)
xlim([0,kMaxMax+1])
grid on
ylabel('$\textrm{Average error}$','interpreter','latex','fontsize',18)
xlabel('$k_{Max}$','interpreter','latex','fontsize',18)
legend('Fourier','Brouwer','fontsize',12)
hold off

figure(3)
clf
errorbar(k,iErrM(2:end),iErrM(2:end)-prctile(iErrMat(2:end),20),prctile(iErrMat(2:end),80)-iErrM(2:end),'o','linewidth',2)
hold on
plot(0:(kMaxMax+1),iErrM(1)*ones(1,length(k)+2),'--','linewidth',2)
xticks(1:kMaxMax)
xlim([0,kMaxMax+1])
grid on
ylabel('$\textrm{Average error}$','interpreter','latex','fontsize',18)
xlabel('$k_{Max}$','interpreter','latex','fontsize',18)
legend('Fourier','Brouwer','fontsize',12)
hold off

figure(4)
clf
errorbar(k,raanErrM(2:end),raanErrM(2:end)-prctile(raanErrMat(2:end),20),prctile(raanErrMat(2:end),80)-raanErrM(2:end),'o','linewidth',2)
hold on
plot(0:(kMaxMax+1),raanErrM(1)*ones(1,length(k)+2),'--','linewidth',2)
xticks(1:kMaxMax)
xlim([0,kMaxMax+1])
grid on
ylabel('$\textrm{Average error}$','interpreter','latex','fontsize',18)
xlabel('$k_{Max}$','interpreter','latex','fontsize',18)
legend('Fourier','Brouwer','fontsize',12)
hold off

figure(5)
clf
errorbar(k,aopErrM(2:end),aopErrM(2:end)-prctile(aopErrMat(2:end),20),prctile(aopErrMat(2:end),80)-aopErrM(2:end),'o','linewidth',2)
hold on
plot(0:(kMaxMax+1),aopErrM(1)*ones(1,length(k)+2),'--','linewidth',2)
xticks(1:kMaxMax)
xlim([0,kMaxMax+1])
grid on
ylabel('$\textrm{Average error}$','interpreter','latex','fontsize',18)
xlabel('$k_{Max}$','interpreter','latex','fontsize',18)
legend('Fourier','Brouwer','fontsize',12)
hold off

figure(6)
clf
errorbar(k,MErrM(2:end),MErrM(2:end)-prctile(MErrMat(2:end),20),prctile(MErrMat(2:end),80)-MErrM(2:end),'o','linewidth',2)
hold on
plot(0:(kMaxMax+1),MErrM(1)*ones(1,length(k)+2),'--','linewidth',2)
xticks(1:kMaxMax)
xlim([0,kMaxMax+1])
grid on
ylabel('$\textrm{Average error}$','interpreter','latex','fontsize',18)
xlabel('$k_{Max}$','interpreter','latex','fontsize',18)
legend('Fourier','Brouwer','fontsize',12)
hold off

%% Plot Relative

aErrRel = aErrMat(:,2:end)-aErrMat(:,1);
eErrRel = eErrMat(:,2:end)-eErrMat(:,1);
iErrRel = iErrMat(:,2:end)-iErrMat(:,1);
raanErrRel = raanErrMat(:,2:end)-raanErrMat(:,1);
aopErrRel = aopErrMat(:,2:end)-aopErrMat(:,1);
MErrRel = MErrMat(:,2:end)-MErrMat(:,1);

f = figure(11);
clf
f.Position(3:4) = [500,250];
errorbar(k,mean(aErrRel),mean(aErrRel)-prctile(aErrRel,10),prctile(aErrRel,90)-mean(aErrRel),'o','linewidth',2)
hold on
plot(0:(kMaxMax+1),zeros(1,length(k)+2),'--k','linewidth',1)
xticks(1:kMaxMax)
xlim([0,kMaxMax+1])
grid on
ylabel('$\frac{\overline{\delta a_F-\delta a_B}}{a\left(0\right)}$','interpreter','latex','fontsize',18)
xlabel('$k_{Max}$','interpreter','latex','fontsize',18)
hold off

f = figure(12);
clf
f.Position(3:4) = [500,250];
errorbar(k,mean(eErrRel),mean(eErrRel)-prctile(eErrRel,10),prctile(eErrRel,90)-mean(eErrRel),'o','linewidth',2)
hold on
plot(0:(kMaxMax+1),zeros(1,length(k)+2),'--k','linewidth',1)
xticks(1:kMaxMax)
xlim([0,kMaxMax+1])
grid on
ylabel('$\frac{\overline{\delta e_F-\delta e_B}}{e\left(0\right)}$','interpreter','latex','fontsize',18)
xlabel('$k_{Max}$','interpreter','latex','fontsize',18)
hold off

f = figure(13);
clf
f.Position(3:4) = [500,250];
errorbar(k,mean(iErrRel),mean(iErrRel)-prctile(iErrRel,10),prctile(iErrRel,90)-mean(iErrRel),'o','linewidth',2)
hold on
plot(0:(kMaxMax+1),zeros(1,length(k)+2),'--k','linewidth',1)
xticks(1:kMaxMax)
xlim([0,kMaxMax+1])
grid on
ylabel('$\overline{\delta i_F-\delta i_B} \left[rad\right]$','interpreter','latex','fontsize',18)
xlabel('$k_{Max}$','interpreter','latex','fontsize',18)
hold off

f = figure(14);
clf
f.Position(3:4) = [500,250];
errorbar(k,mean(raanErrRel),mean(raanErrRel)-prctile(raanErrRel,10),prctile(raanErrRel,90)-mean(raanErrRel),'o','linewidth',2)
hold on
plot(0:(kMaxMax+1),zeros(1,length(k)+2),'--k','linewidth',1)
xticks(1:kMaxMax)
xlim([0,kMaxMax+1])
grid on
ylabel('$\overline{\delta \Omega_F-\delta \Omega_B} \left[rad\right]$','interpreter','latex','fontsize',18)
xlabel('$k_{Max}$','interpreter','latex','fontsize',18)
hold off

f = figure(15);
clf
f.Position(3:4) = [500,250];
errorbar(k,mean(aopErrRel),mean(aopErrRel)-prctile(aopErrRel,10),prctile(aopErrRel,90)-mean(aopErrRel),'o','linewidth',2)
hold on
plot(0:(kMaxMax+1),zeros(1,length(k)+2),'--k','linewidth',1)
xticks(1:kMaxMax)
xlim([0,kMaxMax+1])
grid on
ylabel('$\overline{\delta \omega_F-\delta \omega_B} \left[rad\right]$','interpreter','latex','fontsize',18)
xlabel('$k_{Max}$','interpreter','latex','fontsize',18)
hold off

f = figure(16);
clf
f.Position(3:4) = [500,250];
errorbar(k,mean(MErrRel),mean(MErrRel)-prctile(MErrRel,10),prctile(MErrRel,90)-mean(MErrRel),'o','linewidth',2)
hold on
plot(0:(kMaxMax+1),zeros(1,length(k)+2),'--k','linewidth',1)
xticks(1:kMaxMax)
xlim([0,kMaxMax+1])
grid on
ylabel('$\overline{\delta M_F-\delta M_B} \left[rad\right]$','interpreter','latex','fontsize',18)
xlabel('$k_{Max}$','interpreter','latex','fontsize',18)
hold off