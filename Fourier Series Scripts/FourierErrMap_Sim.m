clear
%% Define Parameters
dataFolder = 'C:\Users\User\Dropbox\Doc Fourier Data\Error Mapping';
dbPath = 'C:\Users\User\Dropbox'; % ASRI
primary = earth();
kMax = 5;
nOrb = 1;
dT = 100; % sec
depritFlag = 1;
nT = 80;
% Region Params
nInc = 360;
nEcc = 100;
nMonte = 1000; % 6000 trials is about 1.3 minute (not parallel)
incRange = linspace(0,90,nInc);
% eccRange = linspace(0.01,0.5,nEcc);
eccRange = logspace(-2,-1,nEcc);
maxSma = 25000;

%% Initialize Error Tensors
errTenF = inf(nEcc,nInc,6);
errTenF2 = inf(nEcc,nInc,6);
errTenB = inf(nEcc,nInc,6);
errTenD = inf(nEcc,nInc,6);

bTime = 0;
cTime = 0;
fTime = 0;
f2Time = 0;
dTime = 0;
%% Loops
disp(['Starting Mapping' newline 'Estimated runtime: ' num2str(nInc*nEcc*nMonte/4500/4/60) 'h'])
totalTime = tic;
parfor iEcc = 1:nEcc
    ecc = eccRange(iEcc);
    minSma = (primary.Re+100)/(1-ecc); % Can change to not go through the atmosphere
    for iInc = 1:nInc
        inc = incRange(iInc);
        % Average Error Vectors
        errVecB = zeros(6,1);
        errVecF = zeros(6,1);
        errVecF2 = zeros(6,1);
        errVecD = zeros(6,1);
        for iTry = 1:nMonte
            % Get OE
            rNum = rand(3,1);
            sma = minSma + rNum(1)*(maxSma-minSma);
            ran = rNum(2)*360;
            aop = rNum(3)*360;
            man = 0;
            oe = [sma,ecc,inc,ran,aop,man];
            
            
            % Define Sat
            Sat = SingleSat(oe,earth());
            Prop = Propagator(Sat);
            % Get Time Vector - prop Deprit first
            T = 2*pi*sqrt(oe(1)^3/Sat.primary.mu);
            if depritFlag
                try
                    tic
                    [t,oeD] = Prop.PropOeDeprit(nT,nOrb);
                    
                    testT = toc;
                    dTime = dTime + testT;
                    oeD = oeD.';
                    if any(diff(t)<0)
                        errorIFTTT(dbPath,'non-monotonic time!');
                    end
                catch
                    t = linspace(0,nOrb*T,nT);
                    oeD = nan(6,80);
                end
            else
                t = 0:dT:nOrb*T; %#ok<UNRCH>
                oeD = nan(6,80);
            end
            % Prop Numerical
            tic
            [~,oeC] = Prop.PropOeOsc(t);
            testT = toc;
            cTime = cTime + testT;
            oeC = oeC.';
            
            try
                % Prop Brouwer
                tic
                [~,OeM] = Prop.PropOeMeanShort(t);
                oeB = me2oscSP(OeM.');
                oeB(6,:) = 180/pi*unwrap(pi/180*oeB(6,:));
                if oeB(6,1) > 180
                    oeB(6,:) = oeB(6,:) - 360;
                end
                testT = toc;
                bTime = bTime + testT;
                errB = abs(oeC-oeB);
                errB = [errB(1,:)/oe(1);errB(2,:)/oe(2);errB(3:end,:)*pi/180];
                errVecB = errVecB + trapz(t.',errB,2)/t(end);
            catch
                % Error remains infinite
                errVecB = errVecB + inf(6,1);
            end
            % Prop Fourier
            tic
            [~,oeF] = Prop.PropOeFourier(t,kMax);
            testT = toc;
            fTime = fTime + testT;
            oeF = oeF.';
            
            % Prop Fourier 2nd order
            tic
            [~,oeF2] = Prop.PropOeFourier2Ord(t,kMax);
            testT = toc;
            f2Time = f2Time + testT;
            oeF2 = oeF2.';
            
            % Fourier Error
            errF = abs(oeC-oeF);
            errF = [errF(1,:)/oe(1);errF(2,:)/oe(2);errF(3:end,:)*pi/180];
            errVecF = errVecF + trapz(t.',errF,2)/t(end);
            
            % Fourier Error 2nd order
            errF2 = abs(oeC-oeF2);
            errF2 = [errF2(1,:)/oe(1);errF2(2,:)/oe(2);errF2(3:end,:)*pi/180];
            errVecF2 = errVecF2 + trapz(t.',errF2,2)/t(end);
            
            % Deprit Error
            errD = abs(oeC-oeD);
            errD = [errD(1,:)/oe(1);errD(2,:)/oe(2);errD(3:end,:)*pi/180];
            errVecD = errVecD + trapz(t.',errD,2)/t(end);
            
        end
        % Average
        errVecB = errVecB/nMonte;
        errVecF = errVecF/nMonte;
        errVecF2 = errVecF2/nMonte;
        errVecD = errVecD/nMonte;
        % Assign errors
        errTenB(iEcc,iInc,:) = errVecB;
        errTenF(iEcc,iInc,:) = errVecF;
        errTenF2(iEcc,iInc,:) = errVecF2;
        errTenD(iEcc,iInc,:) = errVecD;
    end
end
eTime = toc(totalTime);
reportIFTTT(dbPath,eTime);
%% Save Data
MapData = struct();
MapData.eccRange = eccRange;
MapData.incRange = incRange;
MapData.maxSma = maxSma;
MapData.kMax = kMax;
MapData.nOrb = nOrb;
MapData.dT = dT;
MapData.nT = nT;
MapData.nMonte = nMonte;

MapData.eTime = eTime;
MapData.cTime = cTime;
MapData.dTime = dTime;
MapData.fTime = fTime;
MapData.f2Time = f2Time;
MapData.bTime = bTime;


MapData.errTenF = errTenF;
MapData.errTenF2 = errTenF2;
MapData.errTenB = errTenB;
MapData.errTenD = errTenD;

c = clock;
save([dataFolder '\ErrMaps_' num2str(c(3)) '-' num2str(c(2)) '-' ...
    num2str(c(1)) '_' num2str(c(4)) '-' num2str(c(5)), '.mat'],'MapData');
