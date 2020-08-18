paretoFolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Final Results';
conFolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\LatticeDef v1';
latEm = 40;
elevMin = 5;
load([paretoFolder '/ParetoList_Lat_' num2str(latEm) '_LFC Ex_hASplit.mat']);
hAList = [900, 1000];
dDay = 5;


for iHA = 1:2
    pdopArray = cell(length(paretoSet.nPlanes{iHA}),2);
    for iCon = 1:length(paretoSet.nPlanes{iHA})
        load([conFolder '\LatticeExSol_Lat_' num2str(latEm) '_nSats_' ...
            num2str(paretoSet.nSats{iHA}(iCon)) '_hA_' num2str(hAList(iHA)) '.mat']);
        iPar = find(ExSol.archMat(1,:) == paretoSet.nPlanes{iHA}(iCon));
        
        Arch = struct();
        Arch.nSats = paretoSet.nSats{iHA}(iCon);
        Arch.nPlanes = ExSol.archMat(1,iPar);
        Arch.nAops = ExSol.archMat(2,iPar);
        Arch.nSatsPerAop = ExSol.archMat(3,iPar);
        
        Phase = struct();
        Phase.nC1 = ExSol.phaseMat(1,iPar);
        Phase.nC2 = ExSol.phaseMat(2,iPar);
        Phase.nC3 = ExSol.phaseMat(3,iPar);
        
        Orbit = ExSol.orbits{iPar};
        
        InitCon = ExSol.inits{iPar};
        LcLong = LatticeConstellation(Arch,Phase,Orbit,InitCon);
        PropLong = Propagator(LcLong);
        
        aopRate = 3/4*LcLong.J2*sqrt(LcLong.mu)*LcLong.Re^2*LcLong.sma^(-7/2)*...
            (1-LcLong.ecc^2)^-2*(5*cosd(LcLong.inc)^2-1);
        
        aopTime = ceil(2*pi/abs(aopRate)/86164);
        
        dayVec = 0:dDay:aopTime;
        pdopDay = nan(1,length(dayVec));
        longTimeVec = 86164*dayVec;
        
        [~, mOeDay] = PropLong.PropOeMean(longTimeVec);
        
        for iDay = 1:length(dayVec)
            initDay.M1 = InitCon.M1 + mOeDay(6,iDay);
            initDay.raan1 = InitCon.raan1 + mOeDay(4,iDay);
            initDay.aop1 = InitCon.aop1 + mOeDay(5,iDay);
            
            LcDay = LatticeConstellation(Arch,Phase,Orbit,initDay);
            propDay = Propagator(LcDay,1e-6,1e-6);
            [propTime,propState] = propDay.PropEciJ2(0:100:86164);
            
            pdopC = TdoaPdopVec(propState,propTime,latEm,0,0,elevMin);
            pdopN = TdoaPdopVec(propState,propTime,latEm+5,0,0,elevMin);
            pdopS = TdoaPdopVec(propState,propTime,latEm-5,0,0,elevMin);
            pdop = [pdopN.';pdopC.';pdopS.'];
            pdop(pdop>1000) = 1000;
            pdop(isnan(pdop)) = 1000;
            pdopDay(iDay) = mean(trapz(propTime,pdop,2)/propTime(end));
        end
        pdopArray{iCon,iHA} = pdopDay;
    end
end

%%
close all
for iHA = 1:2
    for iCon = 1:length(paretoSet.nPlanes{iHA})
        figure(iHA)
        hold on 
        semilogy((0:length(pdopArray{iCon,iHA})-1)*dDay,pdopArray{iCon,iHA},'-o')
        hold off
    end
end

