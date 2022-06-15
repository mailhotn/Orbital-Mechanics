clear
%% Define Parameters
dataFolder = 'C:\Users\User\Dropbox\Fourier Data\Error Mapping';
primary = earth();
kMax = 4;
nOrb = 1;
dT = 100; % sec
% Region Params
nInc = 50;
nEcc = 20;
nMonte = 1000;
incRange = linspace(60,70,nInc);
eccRange = linspace(0.5,0.7,nEcc);
maxSma = 15000;

%% Initialize Error Tensors
errTenF = inf(nEcc,nInc,6);
errTenB = inf(nEcc,nInc,6);

%% Loops
for iEcc = 1:nEcc
    ecc = eccRange(iEcc);
    minSma = (primary.Re+1)/(1-ecc); % Can change to not go through the atmosphere
    for iInc = 1:nInc
        inc = incRange(iInc);
        % Average Error Vectors
        errVecB = zeros(6,1);
        errVecF = zeros(6,1);
        parfor iTry = 1:nMonte
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
            T = 2*pi*sqrt(oe(1)^3/Sat.primary.mu);
            t = 0:dT:nOrb*T;
            
            % Prop Numerical
            [~,oeC] = Prop.PropOeOsc(t);
            oeC = oeC.';
            try
                % Prop Brouwer
                [~,OeM] = Prop.PropOeMeanFast(t);
                oeB = me2osc(OeM.');
                errB = abs(oeC-oeB);
                errB = [errB(1,:)/oe(1);errB(2,:)/oe(2);errB(3:end,:)*pi/180];
                errVecB = errVecB + trapz(t.',errB,2)/t(end);
            catch ME
                % Error remains infinite
                errVecB = errVecB + inf(6,1);
            end
            % Prop Fourier
            [~,oeF] = Prop.PropOeFourier2(t,kMax);
            
            % Fourier Error          
            errF = abs(oeC-oeF);
            errF = [errF(1,:)/oe(1);errF(2,:)/oe(2);errF(3:end,:)*pi/180];
            errVecF = errVecF + trapz(t.',errF,2)/t(end);
            
        end
        % Average
        errVecB = errVecB/nMonte;
        errVecF = errVecF/nMonte;
        % Assign errors
        errTenB(iEcc,iInc,:) = errVecB;
        errTenF(iEcc,iInc,:) = errVecF;
    end
end

%% Save Data
MapData = struct();
MapData.eccRange = eccRange;
MapData.incRange = incRange;
MapData.maxSma = maxSma;
MapData.kMax = kMax;
MapData.nOrb = nOrb;
MapData.dT = dT;
MapData.nMonte = nMonte;

MapData.errTenF = errTenF;
MapData.errTenB = errTenB;

c = clock;
save([dataFolder '\ErrMaps_' num2str(c(3)) '-' num2str(c(2)) '-' num2str(c(1)) '_' num2str(c(4)) '-' num2str(c(5)), '.mat'],'MapData');
