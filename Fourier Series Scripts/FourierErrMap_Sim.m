clear
%% Define Parameters
dataFolder = 'C:\Users\User\Dropbox\Fourier Data\Error Mapping';
dbPath = 'C:\Users\User\Dropbox'; % ASRI
primary = earth();
kMax = 4;
nOrb = 1;
dT = 100; % sec
depritFlag = 1;
nT = 80;
% Region Params
nInc = 180;
nEcc = 25;
nMonte = 200; % 200 trials is about 1 minute
incRange = linspace(0.1,90,nInc);
eccRange = linspace(0.001,0.5,nEcc);
maxSma = 25000;

%% Initialize Error Tensors
errTenF = inf(nEcc,nInc,6);
errTenB = inf(nEcc,nInc,6);
errTenD = inf(nEcc,nInc,6);

bTime = 0;
cTime = 0;
fTime = 0;
dTime = 0;
%% Loops
tic
for iEcc = 1:nEcc
    ecc = eccRange(iEcc);
    minSma = (primary.Re+100)/(1-ecc); % Can change to not go through the atmosphere
    for iInc = 1:nInc
        inc = incRange(iInc);
        % Average Error Vectors
        errVecB = zeros(6,1);
        errVecF = zeros(6,1);
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
                [~,OeM] = Prop.PropOeMeanFast(t);
                oeB = me2osc(OeM.');
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
            
            % Fourier Error
            errF = abs(oeC-oeF);
            errF = [errF(1,:)/oe(1);errF(2,:)/oe(2);errF(3:end,:)*pi/180];
            errVecF = errVecF + trapz(t.',errF,2)/t(end);
            
            % Deprit Error
            errD = abs(oeC-oeD);
            errD = [errD(1,:)/oe(1);errF(2,:)/oe(2);errF(3:end,:)*pi/180];
            errVecD = errVecD + trapz(t.',errD,2)/t(end);
            
        end
        % Average
        errVecB = errVecB/nMonte;
        errVecF = errVecF/nMonte;
        errVecD = errVecD/nMonte;
        % Assign errors
        errTenB(iEcc,iInc,:) = errVecB;
        errTenF(iEcc,iInc,:) = errVecF;
        errTenD(iEcc,iInc,:) = errVecD;
    end
end
eTime = toc;
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

MapData.errTenF = errTenF;
MapData.errTenB = errTenB;
MapData.errTenD = errTenD;

c = clock;
save([dataFolder '\ErrMaps_' num2str(c(3)) '-' num2str(c(2)) '-' num2str(c(1)) '_' num2str(c(4)) '-' num2str(c(5)), '.mat'],'MapData');
