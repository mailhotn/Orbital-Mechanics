clear
dataFolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Final Results';
nameList = {'Walker Ex','LFC Ex','LFC Ga'};
folderList = {...
    'C:\Users\User\Dropbox\Walker Optimization Data\Previous Optimization Runs\Latticified Walker multilat multiInc collision';...
    'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\LatticeDef v1';...
    'C:\Users\User\Dropbox\Lattice Optimization Data\GA Standard\Previous Runs\Version 6 - definitive';...
    };
latList = 30:10:60;
% latEm = 30;

prct = 95;

WalkerExPhase = [];
LatticeExPhase = [];
LatticeGaPhase = [];

WalkerExOrbit = [];
LatticeExOrbit = [];
LatticeGaOrbit = [];

WalkerMinDist = [];
LatticeExMinDist = [];
LatticeGaMinDist = [];

for iLat = 1:length(latList)
    close all
    latEm = latList(iLat);
    for iScenario = 1:2:3
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
                Arch.nPlanes = ExSol.archMat(1,iPar);
                Arch.nAops = ExSol.archMat(2,iPar);
                Arch.nSatsPerAop = ExSol.archMat(3,iPar);
                
                Phase = struct();
                Phase.nC1 = ExSol.phaseMat(1,iPar);
                Phase.nC2 = ExSol.phaseMat(2,iPar);
                Phase.nC3 = ExSol.phaseMat(3,iPar);
                
                Orbit = ExSol.orbits{iPar};
                
                InitCon = ExSol.inits{iPar};
                
            elseif iScenario == 3
                load([folderList{iScenario} '\LatticeGaSol_Lat_' num2str(latEm)...
                    '_nSats_' num2str(paretoSet.nSats(iCon)) '.mat']);
                iPar = find(GaSol.archMat(1,:) == paretoSet.nPlanes(iCon));
                
                LatticeGaPhase = [LatticeGaPhase, GaSol.phaseMat(1:3,iPar)];
                LatticeGaOrbit = [LatticeGaOrbit, [GaSol.orbits{iPar}.sma;...
                    GaSol.orbits{iPar}.ecc;GaSol.orbits{iPar}.inc;]];
                
                Arch = struct();
                Arch.nSats = paretoSet.nSats(iCon);
                Arch.nPlanes = GaSol.archMat(1,iPar);
                Arch.nAops = GaSol.archMat(2,iPar);
                Arch.nSatsPerAop = GaSol.archMat(3,iPar);
                
                Phase = struct();
                Phase.nC1 = GaSol.phaseMat(1,iPar);
                Phase.nC2 = GaSol.phaseMat(2,iPar);
                Phase.nC3 = GaSol.phaseMat(3,iPar);
                
                Orbit = GaSol.orbits{iPar};
                
                InitCon = GaSol.inits{iPar};
                
            end
            
            
            LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);
            minDist = CalcMinDist(LC);
            
            figure(iScenario*100 + iCon*10 + 1)
            %         subplot(2,1,1)
            %         PlotGroundTrack(LC);
            %         title([nameList{iScenario} ' ' num2str(Arch.nSats) '/' ...
            %             num2str(Arch.nPlanes) ' ' num2str(Orbit.hA)])
            %
            %         subplot(2,1,2)
            PlotPdopMap(LC,latEm,0,5,true,true,prct);
            
            switch iScenario
                case 1
                    WalkerMinDist = [WalkerMinDist, minDist];
                    name = [num2str(prct) 'Lat' num2str(latEm) 'Wa ' num2str(LC.nSats)...
                        '-' num2str(LC.nPlanes)];
                case 2
                    LatticeExMinDist{iHA} = [LatticeExMinDist{iHA}, minDist];
                    name = [num2str(prct) 'Lat' num2str(latEm) 'Ex900 ' num2str(LC.nSats)...
                        '-' num2str(LC.nPlanes)];
                case 3
                    LatticeGaMinDist = [LatticeGaMinDist, minDist];
                    name = [num2str(prct) 'Lat' num2str(latEm) 'Ga ' num2str(LC.nSats)...
                        '-' num2str(LC.nPlanes)];
            end
            
            print([dataFolder '\Coverage Maps\' name],'-depsc','-painters');
        end
    end
end