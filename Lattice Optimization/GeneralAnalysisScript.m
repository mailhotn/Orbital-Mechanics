%% Define Scenarios for Analysis
folderList = {'C:\Users\User\Dropbox\Walker Optimization Data\Previous Optimization Runs\Walker RGT Ex Search delta inc 10',...
    'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\Lattice Version 3',...
    'C:\Users\User\Dropbox\Lattice Optimization Data\Previous Runs\Lattice O-M Optimal Inc',...
    'C:\Users\User\Dropbox\Lattice Optimization Data'};
markerList = {'*','o','x','s'};
nameList = {'Walker','Lattice Basic','Lattice Optimal 1','Lattice Optimal 2'};
nScenarios = numel(folderList);
%% Set Performance Goals
intTarget = 3;
maxTarget = 10;
covTarget = 99.5;

%% Analyze Scenarios
for iScenario = 1:nScenarios
    load([folderList{iScenario} '\OptParams.mat']);
    nSats = minSats:maxSats;
    maxPdop = nan(numel(eccList),numel(nSats));
    intPdop = nan(numel(eccList),numel(nSats));
    coverage = nan(numel(eccList),numel(nSats));
    for iLat = 1:numel(latList)
        for iEcc = 1:numel(eccList)
            for iSats = 1:numel(nS
                ats)
                load([folderList{iScenario} '\LatticeExSol_Lat_'...
                    num2str(latList(iLat)) '_nSats_' num2str(nSats(iSats)) ...
                    '_ecc_' num2str(eccList(iEcc)) '.mat']);
                maxPdop(iEcc,iSats) = ExSol.maxPdop(ExSol.iOpt);
                intPdop(iEcc,iSats) = ExSol.intPdop(ExSol.iOpt);
                coverage(iEcc,iSats) = ExSol.coverage(ExSol.iOpt);
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
        hold off
    end
end