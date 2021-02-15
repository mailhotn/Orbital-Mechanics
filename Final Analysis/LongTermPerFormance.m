clear
paretoFolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Final Results';
conFolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\LatticeDef v1';
latList = 30:10:60;
% latEm = 50;
elevMin = 5;

hAList = [900, 1000];
dDay = 2;

%% Find Nominal Simulation Length
for iLat = 1:length(latList)
    latEm = latList(iLat);
    load([paretoFolder '/ParetoList_Lat_' num2str(latEm) '_LFC Ex_hASplit.mat']);
    aopTime = 0;
    for iHA = 1:1
        for iCon = 1:length(paretoSet.nPlanes{iHA+1})
            load([conFolder '\LatticeExSol_Lat_' num2str(latEm) '_nSats_' ...
                num2str(paretoSet.nSats{iHA+1}(iCon)) '_hA_' num2str(hAList(iHA)) '.mat']);
            iPar = find(ExSol.archMat(1,:) == paretoSet.nPlanes{iHA+1}(iCon));
            Orbit = ExSol.orbits{iPar};
            primary = earth();
            
            aopRate = 3/4*primary.J2*sqrt(primary.mu)*primary.Re^2*Orbit.sma^(-7/2)*...
                (1-Orbit.ecc^2)^-2*(5*cosd(Orbit.inc)^2-1);
            aopTime = max([aopTime, ceil(2*pi/abs(aopRate)/86164)]);
        end
    end
    aopTime = min([aopTime, 2*365]);
    dayVec = 0:dDay:aopTime;
    longTimeVec = 86164*dayVec;
    
    pdop900 = nan(length(paretoSet.nPlanes{2}),length(dayVec));
    pdop1000 = nan(length(paretoSet.nPlanes{3}),length(dayVec));
    
    nRT900 = nan(1,length(paretoSet.nPlanes{2}));
    nRT1000 = nan(1,length(paretoSet.nPlanes{3}));
    
    leg900 = cell(1);
    leg1000 = cell(1);
    %% Run Stuff
    
    for iHA = 1:1
        
        for iCon = 1:length(paretoSet.nPlanes{iHA+1})
            load([conFolder '\LatticeExSol_Lat_' num2str(latEm) '_nSats_' ...
                num2str(paretoSet.nSats{iHA+1}(iCon)) '_hA_' num2str(hAList(iHA)) '.mat']);
            iPar = find(ExSol.archMat(1,:) == paretoSet.nPlanes{iHA+1}(iCon));
            
            Arch = struct();
            Arch.nSats = paretoSet.nSats{iHA+1}(iCon);
            Arch.nPlanes = ExSol.archMat(1,iPar);
            Arch.nAops = ExSol.archMat(2,iPar);
            Arch.nSatsPerAop = ExSol.archMat(3,iPar);
            
            Phase = struct();
            Phase.nC1 = ExSol.phaseMat(1,iPar);
            Phase.nC2 = ExSol.phaseMat(2,iPar);
            Phase.nC3 = ExSol.phaseMat(3,iPar);
            nSR = gcd(gcd(abs(gcd(Arch.nPlanes,Phase.nC3)*Arch.nSatsPerAop),...
                abs(gcd(Arch.nPlanes,Phase.nC3)*1)),abs(Phase.nC1*1-Arch.nSatsPerAop*14));
            nRT = Arch.nSats/nSR;
            Orbit = ExSol.orbits{iPar};
            
            InitCon = ExSol.inits{iPar};
            LcLong = LatticeConstellation(Arch,Phase,Orbit,InitCon);
            PropLong = Propagator(LcLong);
            
            [~, mOeDay] = PropLong.PropOeMean(longTimeVec);
            
            pdopDay = nan(1,length(dayVec));
            parfor iDay = 1:length(dayVec)
                initDay = struct();
                initDay.M1 = mOeDay(iDay,6);
                initDay.raan1 = mOeDay(iDay,4);
                initDay.aop1 = mOeDay(iDay,5);
                
                LcDay = LatticeConstellation(Arch,Phase,Orbit,initDay);
                propDay = Propagator(LcDay,1e-6,1e-6);
                [propTime,propState] = propDay.PropEciJ2(0:100:86164);
                
                pdop = TdoaPdopVec(propState,propTime,latEm,0,0,elevMin);
                %             pdopN = TdoaPdopVec(propState,propTime,latEm+5,0,0,elevMin);
                %             pdopS = TdoaPdopVec(propState,propTime,latEm-5,0,0,elevMin);
                %             pdop = [pdopN.';pdopC.';pdopS.'];
                pdop(pdop>100) = 100;
                pdop(isnan(pdop)) = 100;
                %             pdopDay(iDay) = trapz(propTime,pdop)/propTime(end);
                pdopDay(iDay) = prctile(pdop,95);
            end
            switch iHA
                case 1
                    pdop900(iCon,:) = pdopDay;
                    nRT900(iCon) = nRT;
                    leg900{iCon} = ['LFC ' num2str(Arch.nSats) '/' num2str(Arch.nPlanes)];
                case 2
                    pdop1000(iCon,:) = pdopDay;
                    nRT1000(iCon) = nRT;
                    leg1000{iCon} = ['LFC ' num2str(Arch.nSats) '/' num2str(Arch.nPlanes)];
            end
        end
    end
    
    %% Plot Stuff
    figure(1)
    semilogy(dayVec,pdop900,'--o')
    xlim([0,dayVec(end)])
    ylim([1,100])
    xlabel('$Day$','interpreter','latex','fontsize',14)
    ylabel('$p_{95}\left(PDOP\right)$','interpreter','latex','fontsize',14)
    legend(leg900,'location','best')
    grid on
    
    print([paretoFolder '\Coverage Maps\Long Term hA900 Lat' num2str(latEm)],...
        '-depsc','-painters');
    
    % figure(2)
    % semilogy(dayVec,pdop1000,'-o')
    % xlim([0,dayVec(end)])
    % ylim([1,100])
    % xlabel('$Day$','interpreter','latex','fontsize',14)
    % ylabel('$p_{95}\left(PDOP\right)$','interpreter','latex','fontsize',14)
    % legend(leg1000,'location','best')
    % grid on
end