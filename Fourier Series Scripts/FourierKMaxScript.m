clear
nSats = 100; % different orbits
kMaxMax = 15;
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
        0.1+0.1*r(2);
        63.4+0.2*(2*r(3)-1);
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

%% Plot
k = 1:kMaxMax;
aStd = std(aErrMat);
eStd = std(eErrMat);
iStd = std(iErrMat);
raanStd = std(raanErrMat);
aopStd = std(aopErrMat);
MStd = std(MErrMat);

aErrM = mean(aErrMat);
eErrM = mean(eErrMat);
iErrM = mean(iErrMat);
raanErrM = mean(raanErrMat);
aopErrM = mean(aopErrMat);
MErrM = mean(MErrMat);

figure(1)
clf
errorbar(k,aErrM(2:end),aStd(2:end),'o')
hold on
plot(k,aErrM(1)*ones(1,length(k)),'--')
hold off

figure(2)
clf
errorbar(k,eErrM(2:end),eStd(2:end),'o')
hold on
plot(k,eErrM(1)*ones(1,length(k)),'--')
hold off

figure(3)
clf
errorbar(k,iErrM(2:end),iStd(2:end),'o')
hold on
plot(k,iErrM(1)*ones(1,length(k)),'--')
hold off

figure(4)
clf
errorbar(k,raanErrM(2:end),raanStd(2:end),'o')
hold on
plot(k,raanErrM(1)*ones(1,length(k)),'--')
hold off

figure(5)
clf
errorbar(k,aopErrM(2:end),aopStd(2:end),'o')
hold on
plot(k,aopErrM(1)*ones(1,length(k)),'--')
hold off

figure(6)
clf
errorbar(k,MErrM(2:end),MStd(2:end),'o')
hold on
plot(k,MErrM(1)*ones(1,length(k)),'--')
hold off

