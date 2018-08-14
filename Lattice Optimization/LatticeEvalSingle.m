
datafolder = 'C:\Users\User\Dropbox\Lattice Optimization Data';
%% Choose Constellation
latGs = 40;
ecc = 0.05;
nSats = 51;
nPlanes = 51;
nC2 = 1;
Arch.nDays = 1;
Arch.nRepeats = 14; 
load ([datafolder '\OptParams.mat']);
load([datafolder '\LatticeExSol_Lat_' num2str(latGs) ...
    '_nSats_' num2str(nSats) '_ecc_' num2str(ecc) '.mat']);
%% Create Constellation
Arch.nSats = nSats;
Arch.nPlanes = nPlanes;
Phase.nC3 = Arch.nPlanes; % Constraint! can be different, would increase number of relative orbits.
Arch.nSatsPerAop = Arch.nDays; % Constraint!  can be different if nDays > 1
Arch.nAops = Arch.nSats/Arch.nPlanes/Arch.nSatsPerAop;
gcdOrbits = gcd(Arch.nPlanes,Phase.nC3);
if gcdOrbits > Arch.nRepeats
    Phase.nC1 = Arch.nRepeats;
else
    l = floor(Arch.nRepeats/gcdOrbits);
    if l*gcdOrbits == Arch.nRepeats
        l = l-1;
    end
    Phase.nC1 = Arch.nRepeats - l*gcdOrbits;
end
Phase.nC2 = nC2;

LC = LatticeConstellation(Arch,Phase,ExSol.Orbit,ExSol.InitCon);
%% Propagate
Prop = Propagator(LC,PropParams.relTol,PropParams.absTol);
[propTime, propState] = Prop.PropEciJ2(PropParams.timeVec);
% Evaluate PDOP
[pdop, satsIs] = TdoaPdopVec(propState,propTime,latGs,0,0,PropParams.elevMin);
%% Plot Stuff
figure(1)
PlotGroundTrack(propState,propTime,0)

figure(2)
yyaxis left
plot(propTime/3600,pdop)
ylabel('PDOP')
yyaxis right
plot(propTime/3600,satsIs)
ylabel('Sats In Sight')
xlabel('Time [hr]')
grid