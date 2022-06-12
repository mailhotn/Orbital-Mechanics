%% Define Parameters
primary = earth();
kMax = 4;
nOrb = 1;
dT = 100; % sec
% Region Params
nInc = 180;
nEcc = 100;
nMonte = 1000;
incRange = linspace(0.4,90,nInc);
eccRange = linspace(0.001,0.3,nEcc);
maxSma = 15000;

% Error Tensors
errTenF = nan(nEcc,nInc,6);
errTenB = nan(nEcc,nInc,6);

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
            % Prop Brouwer
            [~,OeM] = Prop.PropOeMeanFast(t);
            oeB = me2osc(OeM.');
            % Prop Fourier
            [~,oeF] = Prop.PropOeFourier2(t,kMax);
            
            % Errors
            errB = abs(oeC-oeB);
            errB = [errB(1,:)/oe(1);errB(2,:)/oe(2);errB(3:end,:)*pi/180];
            errVecB = errVecB + trapz(t.',errB,2)/t(end);
            
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

%% Plot
[incMesh, eccMesh] = meshgrid(incRange,eccRange);
levelVec = 0.01*[-1:0.01:0,0:0.01,1];

figure(1)
contourf(incMesh,eccMesh,errTenF(:,:,1)-errTenB(:,:,1))
colorbar

figure(2)
contourf(incMesh,eccMesh,errTenF(:,:,2)-errTenB(:,:,2))
colorbar

figure(3)
contourf(incMesh,eccMesh,errTenF(:,:,3)-errTenB(:,:,3))
colorbar

figure(4)
contourf(incMesh,eccMesh,errTenF(:,:,4)-errTenB(:,:,4))
colorbar

figure(5)
contourf(incMesh,eccMesh,errTenF(:,:,5)-errTenB(:,:,5))
colorbar

figure(6)
contourf(incMesh,eccMesh,errTenF(:,:,6)-errTenB(:,:,6))
colorbar