dataFolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Final Results';
% load([datafolder '\OptParams.mat']);
incList = 0:90;
hAList = [0,900,1000];
nRepeats = 14;
nDays = 1;
%% Drag Parameters
sDivM = 1;
cD = 1;
kD = sDivM*cD;
%% Fuel Consumption For All ecc & Inc
delV = nan(length(incList),length(hAList));
thPulse = nan(length(incList),length(hAList));
eccs = nan(length(incList),length(hAList));
smas = nan(length(incList),length(hAList));
for iInc = 1:length(incList)
    for iHA = 1:length(hAList)
        
        hA = hAList(iHA);
        inc = incList(iInc);
        
        [sma, ecc] = CalcRgtSmaApoHeight(inc,hA,nRepeats,nDays);
        eccs(iInc,iHA) = ecc;
        smas(iInc,iHA) = sma;
        
        [rhoP, ~, ~, hScale] = EarthAtmosphere(sma*(1-ecc));
        z = sma*ecc/hScale;
        % Calculations in meters!!! (rho is kg/m^3)
        delSmaF = -2*pi*kD*rhoP*sma^2*1000^2*exp(-z)*...
            ((1 + 3/4*ecc^2 + 21/64*ecc^4)*besseli(0,z) + ...
            (2*ecc + 3/4*ecc^3)*besseli(1,z) + ...
            (3/4*ecc^2 + 7/16*ecc^4)*besseli(2,z) + ...
            1/4*ecc^3*besseli(3,z) + ...
            7/64*ecc^4*besseli(4,z));
        
        delEcc = -2*pi*kD*rhoP*sma*1000*(1-ecc^2)*exp(-z)*...
            ((1/2*ecc + 3/16*ecc^3)*besseli(0,z) + ...
            (1 + 3/8*ecc^2 + 15/64*ecc^4)*besseli(1,z) + ...
            (1/2*ecc + 1/4*ecc^3)*besseli(2,z) + ...
            (1/8*ecc^2 + 15/128*ecc^4)*besseli(3,z) + ...
            1/16*ecc^3*besseli(4,z) + ...
            3/128*ecc^4*besseli(5,z));
        options = optimset('Display','none');
        
        [th, delVF] = fminbnd(@(x)FuelCostElliptical(x,sma,ecc,-delSmaF,-delEcc)...
            ,0,pi,options);
        delV(iInc,iHA) = delVF;
        thPulse(iInc,iHA) = 180/pi*th;
    end
end
delV = delV*nRepeats*365;
figure(1)
plot(incList.',delV(:,1),'-',...
    incList.',delV(:,2),'-.',...
    incList.',delV(:,3),'--','linewidth',2)
grid on
legend('Circular, h_a = h_c','Elliptical, h_a = 900 km','Elliptical, h_a = 1000 km')
xlabel('$\rm{Orbital \: inclination}$','fontsize',14,'interpreter','latex')
xticks(0:10:90)
xticklabels({'0°','10°','20°','30°','40°','50°','60°','70°','80°','90°'})
ylabel('$\mathrm{Annual} \: \Delta v \mathrm{\left[\frac{m}{s}\right]}$','fontsize',14,'interpreter','latex')
xlim([0,90])
ylim([0,25])

print([dataFolder '\Annual fuel usage'],'-depsc','-painters');
% figure(2)
% plot(incList,thPulse)
%% New Version 08/2020
% syms a e th mu da de real
% p = a*(1-e^2);
% r = p/(1+e*cos(th));
% v = sqrt(mu/p)*sqrt(1+e^2+2*e*cos(th));
% dvT = -da*r*v/2/a/(2*a-r);
% dvN = (-de-2/v*(e+cos(th))*dvT)*a*v/r/sin(th);
% 
% dv = simplify(sqrt(dvT^2 + dvN^2));
% dvdth = diff(dv,th);

%% Function Definition
function [J] = FuelCostElliptical(x,sma,ecc,delSma,delEcc)
primary = earth();
mu = primary.mu;
th = x;
p = sma*(1-ecc^2);
r = p./(1+ecc.*cos(th));
v = 1000*sqrt(mu./p)*sqrt(1 + ecc^2 + 2*ecc*cos(th));
sma = 1000*sma;
r = 1000*r;

dVTan = delSma.*r.*v./(2*sma.*(2*sma-r));

dVRad = delEcc*sma*v/r/sin(th) + 2*sma*(ecc+cos(th))/r/sin(th)*dVTan;

J = sqrt(dVTan.^2 + dVRad.^2);
end
