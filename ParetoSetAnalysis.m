clear
dataFolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Final Results';
nameList = {'Walker Ex','Lattice Ex','Lattice Ga'};
folderList = {...
    'C:\Users\User\Dropbox\Walker Optimization Data\Previous Optimization Runs\Latticified Walker multilat multiInc collision';...
    'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\LatticeDef v1';...
    'C:\Users\User\Dropbox\Lattice Optimization Data\GA Standard\Previous Runs\Version 6 - definitive';...
    };
latEm = 30;
WalkerExPhase = [];
LatticeExPhase = [];
LatticeGaPhase = [];
WalkerExOrbit = [];
LatticeExOrbit = [];
LatticeGaOrbit = [];

for iScenario = 1:2
    load([dataFolder '/ParetoList_Lat_' num2str(latEm) '_' nameList{iScenario} '.mat']);
    for iCon = 1:length(paretoSet.nPlanes)
        if iScenario <= 2
            load([folderList{iScenario} '\LatticeExSol_Lat_'...
                num2str(latEm) '_nSats_' num2str(paretoSet.nSats(iCon)) ...
                '_hA_' num2str(paretoSet.hAs(iCon)) '.mat']);
            
            iPar = find(ExSol.archMat(1,:) == paretoSet.nPlanes(iCon));
            if iScenario == 1
                WalkerExPhase = [WalkerExPhase, ExSol.phaseMat(1:3,iPar)];
                WalkerExOrbit = [WalkerExOrbit, [ExSol.orbits{iPar}.sma;...
                    ExSol.orbits{iPar}.ecc;ExSol.orbits{iPar}.inc;]];
            elseif iScenario == 2
                LatticeExPhase = [LatticeExPhase, ExSol.phaseMat(1:3,iPar)];
                LatticeExOrbit = [LatticeExOrbit, [ExSol.orbits{iPar}.sma;...
                    ExSol.orbits{iPar}.ecc;ExSol.orbits{iPar}.inc;]];
            end
            Arch = struct();
            Arch.nSats = paretoSet.nSats(iCon);
            Arch.nPlanes = paretoSet.nPlanes(iCon);
            Arch.nSatsPerAop = 1;
            Arch.nAops = Arch.nSats/Arch.nPlanes/Arch.nSatsPerAop;
            
            Phase = struct();
            Phase.nC1 = ExSol.phaseMat(1,iPar);
            Phase.nC2 = ExSol.phaseMat(2,iPar);
            Phase.nC3 = ExSol.phaseMat(3,iPar);
            
            Orbit = ExSol.orbits{iPar};
            
            InitCon = ExSol.inits{iPar};
            
            LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);
            
            figure(iScenario*100 + iCon*10 + 1)
            subplot(2,1,1)
            PlotGroundTrack(LC);
                        title([nameList{iScenario} ' ' num2str(Arch.nSats) '/' ...
                num2str(Arch.nPlanes) ' ' num2str(Orbit.hA)])
            
            subplot(2,1,2)
            PlotPdopMap(LC,latEm);

            
            
            
        end
    end
    
end
