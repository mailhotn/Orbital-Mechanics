T = 30;
latGs = 80;
datafolder = 'C:\Users\User\Dropbox\Walker Optimization Data';

load([datafolder '\WalkerRgtExSol_Lat_' num2str(latGs)...
            '_T_' num2str(T) '.mat']);
WC = WalkerConstellation(ExSol.nSats,ExSol.optP,ExSol.optF,ExSol.inc,...
    ExSol.alt,ExSol.raan0);
Prop = Propagator(WC,1e-8,1e-9);

[propTime, propState] = Prop.PropEciJ2(ExSol.PropParams.timeVec);
[pdop, satsIs] = TdoaPdopVec(propState,propTime,ExSol.latGs,0,0, ...
    ExSol.PropParams.elevMin);
oeM = WC.InitialOeMean;

figure(1)
yyaxis left
plot(propTime/3600,pdop)
ylabel('PDOP')
yyaxis right
plot(propTime/3600,satsIs)
ylabel('Sats In Sight')
xlabel('Time [hr]')
grid

figure(2)
plot(oeM(4,:),oeM(6,:),'o')
xlabel('\Omega_0')
ylabel('M_0')
axis equal
xlim([0,360])
ylim([0,360])