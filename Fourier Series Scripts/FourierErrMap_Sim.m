clear
TurnOffPCWhenDone = true;
%% Define Parameters
dataFolder = 'C:\Users\User\Google Drive\Doc Data\Error Mapping';
dbPath = 'C:\Users\User\Google Drive'; % ASRI
primary = earth();
k1 = 5;
k2 = 5;
nOrb = 1;
dT = 100; % sec only used if Deprit is not being tested
depritFlag = 0; % Probably forget about this part?
nT = 80; % Only used if Deprit is being tested
% Region Params
nInc = 361;
nEcc = 100;
nMonte = 800; % 10000 trials is about 3.63 minute (not parallel)

% Region parameters for speed test - 6000 runs
% nInc = 12;
% nEcc = 10;
% nMonte = 200;

incRange = linspace(0,90,nInc);
% eccRange = linspace(0.01,0.5,nEcc);
eccRange = logspace(log10(0.002),log10(0.1),nEcc);
maxSma = 10000;

%% Initialize Error Tensors
errTenF = inf(nEcc,nInc,6);
errTenF2 = inf(nEcc,nInc,6);
errTenB = inf(nEcc,nInc,6);
errTenD = inf(nEcc,nInc,6);
% cartesian errors
errTenRswF = inf(nEcc,nInc,6);
errTenRswB = inf(nEcc,nInc,6);

bTime = 0;
cTime = 0;
fTime = 0;
f2Time = 0;
dTime = 0;
%% Loops
disp(datetime('now'))
disp(['Starting Mapping' newline 'Estimated runtime: ' num2str(nInc*nEcc*nMonte*3.64/10000/4/60) 'h'])

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
        errVecRswF = zeros(6,1);
        errVecRswB = zeros(6,1);
        for iTry = 1:nMonte
            % Get OE
            rNum = rand(4,1);
            sma = minSma + rNum(1)*(maxSma-minSma);
            ran = rNum(2)*360;
            aop = rNum(3)*360;
            man = rNum(4)*360;
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
                t = 0:dT:nOrb*T;
                oeD = nan(6,80);
            end
            % Prop Numerical
            tic
            [~,oeC] = Prop.PropOeOsc3(t);
            testT = toc;
            cTime = cTime + testT;
            oeC = oeC.';
            eciC = oe2eci(oeC,primary,'me');

            try
                % Prop Brouwer
                tic
                [~,OeM] = Prop.PropOeMeanShort(t);
                oeB = me2oscSP(OeM.');
                oeB(6,:) = 180/pi*unwrap(pi/180*oeB(6,:));
                eciB = oe2eci(oeB,primary,'me');
                % if oeB(6,1) > 180
                %     oeB(6,:) = oeB(6,:) - 360;
                % end
                testT = toc;
                bTime = bTime + testT;
                errB = abs(oeC-oeB);
                errB = [errB(1,:)/oe(1);errB(2,:)/oe(2);errB(3:end,:)*pi/180];
                errVecB = errVecB + trapz(t.',errB,2)/t(end);

                % Cartesian
                errRswB = abs(eci2rsw(eciB-eciC,oeC));
                errVecRswB = errVecRswB + trapz(t.',errRswB,2)/t(end);
            catch
                % Error remains infinite
                errVecB = errVecB + inf(6,1);
            end

            try
                % Prop Fourier - k1
                tic
                [~,oeF] = Prop.PropOeFourier2Ord(t,k1);
                eciF = oe2eci(oeF,primary,'me');
                testT = toc;
                fTime = fTime + testT;
                oeF = oeF.';

                % Fourier Error
                errF = abs(oeC-oeF);
                errF = [errF(1,:)/oe(1);errF(2,:)/oe(2);errF(3:end,:)*pi/180];
                errVecF = errVecF + trapz(t.',errF,2)/t(end);

                % Cartesian
                errRswF = abs(eci2rsw(eciF-eciC,oeC));
                errVecRswF = errVecRswF + trapz(t.',errRswF,2)/t(end);
            catch
                errVecF = errVecF+inf(6,1);
            end

            if k1~=k2
                try
                    % Prop Fourier - k2
                    tic
                    [~,oeF2] = Prop.PropOeFourier2Ord(t,k2);
                    testT = toc;
                    f2Time = f2Time + testT;
                    oeF2 = oeF2.';
                    % Fourier Error k2

                    errF2 = abs(oeC-oeF2);
                    errF2 = [errF2(1,:)/oe(1);errF2(2,:)/oe(2);errF2(3:end,:)*pi/180];
                    errVecF2 = errVecF2 + trapz(t.',errF2,2)/t(end);
                catch
                    errVecF2 = errVecF2+inf(6,1);
                end
            end

            if depritFlag
                % Deprit Error
                errD = abs(oeC-oeD); %#ok<UNRCH>
                errD = [errD(1,:)/oe(1);errD(2,:)/oe(2);errD(3:end,:)*pi/180];
                errVecD = errVecD + trapz(t.',errD,2)/t(end);
            end

        end
        % Average
        errVecB = errVecB/nMonte;
        errVecF = errVecF/nMonte;
        errVecF2 = errVecF2/nMonte;
        errVecD = errVecD/nMonte;
        errVecRswF = errVecRswF/nMonte;
        errVecRswB = errVecRswB/nMonte;
        % Assign errors
        errTenB(iEcc,iInc,:) = errVecB;
        errTenF(iEcc,iInc,:) = errVecF;
        errTenF2(iEcc,iInc,:) = errVecF2;
        errTenD(iEcc,iInc,:) = errVecD;
        errTenRswF(iEcc,iInc,:) = errVecRswF;
        errTenRswB(iEcc,iInc,:) = errVecRswB;
    end
end
eTime = toc(totalTime);
reportIFTTT(dbPath,eTime);
%% Save Data
MapData = struct();
MapData.eccRange = eccRange;
MapData.incRange = incRange;
MapData.maxSma = maxSma;
MapData.k1 = k1;
MapData.k2 = k2;
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
MapData.errTenRswF = errTenRswF;
MapData.errTenRswB = errTenRswB;

c = clock;
save([dataFolder '\ErrMaps_' num2str(c(3)) '-' num2str(c(2)) '-' ...
    num2str(c(1)) '_' num2str(c(4)) '-' num2str(c(5)), '.mat'],'MapData');
%% PC Shutdown!
if TurnOffPCWhenDone
    system('shutdown /s');
end
