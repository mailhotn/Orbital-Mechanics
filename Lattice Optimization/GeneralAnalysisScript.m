%% Define Scenarios for Analysis
folderList = {'C:\Users\User\Dropbox\Walker Optimization Data\Previous Optimization Runs\Latticified Walker',...
    'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\Apogee Height, delta inc',...
    'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\Apogee Height, opt inc(pre-fix)',...
    'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\Apogee Height, opt inc'};
markerList = {'*','o','x','s'};
nameList = {'Walker','Lattice Basic','Lattice Optimal 1','Lattice Optimal 2'};
nScenarios = numel(folderList);
%% Set Performance Goals
intTarget = 10;
maxTarget = 1000;
covTarget = 99;
p90Target = 3;

%% Analyze Scenarios
nLeg = 0;
for iScenario = 1:nScenarios
    load([folderList{iScenario} '\OptParams.mat']);
    hAList = 0;
    nSats = minSats:maxSats;
    maxPdop = nan(numel(hAList),numel(nSats));
    intPdop = nan(numel(hAList),numel(nSats));
    coverage = nan(numel(hAList),numel(nSats));
    p90 = nan(numel(hAList),numel(nSats));
    for iLat = 1:numel(latList)
        for iHA = 1:numel(hAList)
            for iSats = 1:numel(nSats)
                load([folderList{iScenario} '\LatticeExSol_Lat_'...
                    num2str(latList(iLat)) '_nSats_' num2str(nSats(iSats)) ...
                    '_hA_' num2str(hAList(iHA)) '.mat']);
                maxPdop(iHA,iSats) = ExSol.maxPdop(ExSol.iOpt);
                intPdop(iHA,iSats) = ExSol.intPdop(ExSol.iOpt);
                coverage(iHA,iSats) = ExSol.coverage(ExSol.iOpt);
                p90(iHA,iSats) = ExSol.p90(ExSol.iOpt);
            end
        end
        % Plot Results for Latitude
        figure(iLat*10 + 1)
        hold on
        semilogy(nSats,intPdop,markerList{iScenario})
        title(['Integral of PDOP for \phi_0 = ' num2str(latList(iLat))])
        xlabel('# Sats')
        ylabel('$\frac{1}{T}\int^{T}_{0}{PDOPdt}$','interpreter','latex','fontsize',12)
        grid on
        set(gca, 'YScale', 'log')
        leg = get(gca,'Legend');
        if isempty(leg)
            leg = legend(nameList{iScenario});
        else
            leg.String(nLeg+1:nLeg + numel(hAList)) = nameList(iScenario);
        end
        hold off
    end
    nLeg = nLeg + numel(hAList);
end