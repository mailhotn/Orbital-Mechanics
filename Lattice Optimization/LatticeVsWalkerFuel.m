datafolder = ['C:\Users\User\Dropbox\Lattice Optimization Data']; %#ok<NBRAK>
load([datafolder '\OptParams.mat']);

%% Drag Parameters
sDivM = 1;
cD = 1;
kD = sDivM*cD;
%% Fuel Consumption For All ecc & Inc
delV = nan(length(latList),length(eccList));
for iLat = 1:length(latList)
    for iEcc = 1:length(eccList)
        ecc = eccList(iEcc);
        latGs = latList(iLat);
        inc = min([90,latList(iLat)+delInc]);
        sma = CalcRgtSma(ecc,inc,nRepeats,nDays);
        [rhoP, ~, ~, hScale] = EarthAtmosphere(sma*(1-ecc));
        z = sma*ecc/hScale;
        delSmaF = -2*pi*kD*rhoP*sma^2*exp(-z)*...
            ((1 + 3/4*ecc^2 + 21/64*ecc^4)*besseli(0,z) + ...
            (2*ecc + 3/4*ecc^3)*besseli(1,z) + ...
            (3/4*ecc^2 + 7/16*ecc^4)*besseli(2,z) + ...
            1/4*ecc^3*besseli(3,z) + ...
            7/64*ecc^4*besseli(4,z));
        
        delEcc = -2*pi*kD*rhoP*sma*(1-ecc^2)*exp(-z)*...
            ((1/2*ecc + 3/16*ecc^3)*besseli(0,z) + ...
            (1 + 3/8*ecc^2 + 15/64*ecc^4)*besseli(1,z) + ...
            (1/2*ecc + 1/4*ecc^3)*besseli(2,z) + ...
            (1/8*ecc^2 + 15/128*ecc^4)*besseli(3,z) + ...
            1/16*ecc^3*besseli(4,z) + ...
            3/128*ecc^4*besseli(5,z));
        options = optimset('Display','none');
        
        [~, delVF] = fminbnd(@(x)FuelCostElliptical(x,sma,ecc,-delSmaF,-delEcc)...
            ,0,2*pi,options);
        delV(iLat,iEcc) = delVF;
    end
end

%% Function Definition
function [J] = FuelCostElliptical(x,sma,ecc,delSma,delEcc)
primary = earth();
mu = primary.mu;
th = x;
p = sma*(1-ecc^2);
r = p./(1+ecc.*cos(th));
vTheta = sqrt(mu./p).*(1 + ecc.*cos(th));
vRadial = sqrt(mu./p).*ecc.*sin(th);
vMag = sqrt(vTheta.^2 + vRadial.^2);
dVTan = delSma.*r.*vMag./(2*sma.*(2*sma-r));
dVRad = sma.*vMag./r.*(delEcc - 2./vMag.*(ecc + cos(th)).*dVTan);
J = sqrt(dVTan.^2 + dVRad.^2);
end