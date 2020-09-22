%% Generate example of too regional coverage
clear
dataFolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\LatticeDef v1';
load([dataFolder '\LatticeExSol_Lat_40_nSats_61_hA_900.mat']);
latEm = 40;
Arch.nSats = 61;
Arch.nPlanes = 61;
Arch.nAops = 1;
Arch.nSatsPerAop = 1;

Orbit = ExSol.orbits{end};

Phase.nC1 = ExSol.phaseMat(1,end);
Phase.nC2 = ExSol.phaseMat(2,end);
Phase.nC3 = ExSol.phaseMat(3,end);

InitCon = ExSol.inits{end};

LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);
time = 0:1:86164;
Prop = Propagator(LC);
[~, xEci] = Prop.PropEciJ2(time);

primary = earth();
gmst = primary.we*time;
nSats = size(xEci,2)/6;
tracks = nan(2*nSats,length(gmst));
for iTime = 1:length(gmst)
    xEcef = eci2ecef(reshape(xEci(iTime,:).',6,nSats),gmst(iTime));
    rEcef = xEcef(1:3,:)./sqrt(dot(xEcef(1:3,:),xEcef(1:3,:),1));
    lla   = ecef2lla(rEcef.'*primary.Re,0,primary.Re).';
    tracks(:,iTime) = reshape(lla(1:2,:),2*nSats,1);
end

figure(1)
PlotPdopMap(LC,latEm,0,5,true,true);
hold on
for iSat = 1:1
        plot(wrapTo180(tracks(2*iSat,:)),tracks(2*iSat-1,:),'.w');
end
hold off
print([dataFolder '\TooRegional'],'-depsc','-painters');
    

