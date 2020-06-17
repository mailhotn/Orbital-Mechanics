% Useless before restoring inc data
clear
datafolder =     'C:\Users\User\Dropbox\Lattice Optimization Data\GA Standard\Previous Runs\Version 4 - dT 100';
load([datafolder '\OptParams']);
%% Set Performance Goals
intTarget = 1e10;
maxTarget = 1e10;
covTarget = 95;
p90Target = 5;

%%
optdInc = nan(1,length(latList));
stdInc = nan(1,length(latList));
for iLat = 1:numel(latList)
    incList = [];
    hAList = [];
%     haList = nan(1,maxSats-minSats+1);
    for nSats = minSats:maxSats
        load([datafolder '\LatticeGaSol_Lat_' num2str(latList(iLat)) '_nSats_' ...
            num2str(nSats) '.mat']);
        for iCon = 1:size(GaSol.archMat,2)
            achieveInt = mean(GaSol.intPdop,1) < intTarget;
            achieveMax = mean(GaSol.maxPdop,1) < maxTarget;
            achieveCov = mean(GaSol.coverage,1) > covTarget;
            achievep90 = mean(GaSol.p90,1) < p90Target;
            if achieveInt(iCon) & achieveMax(iCon) & achieveCov(iCon) & achievep90(iCon)
                incList = [incList, GaSol.orbits{iCon}.inc];
                hAList = [hAList, GaSol.orbits{iCon}.hA];
            end
        end
    end
    figure(iLat*10+1)
    hist(incList-latList(iLat),12)
    title(['\phi_0 = ' num2str(latList(iLat)) '°'])
    xlabel('i - \phi_0 [°]')
    ylabel('# Constellations')
    optdInc(iLat) = mean(incList)-latList(iLat);
    stdInc(iLat) = std(incList);
%     figure(iLat*10+2)
%     plot(1:sum(hAList==0),incList(hAList==0)-latList(iLat),'o',...
%         1:sum(hAList==900),incList(hAList==900)-latList(iLat),'o',...
%         1:sum(hAList==1000),incList(hAList==1000)-latList(iLat),'o');
%     figure(iLat*10+2)
%     plot(minSats:maxSats,haList,'o')
%     mean(incList-latList(iLat))
%     std(incList-latList(iLat))
end
save([datafolder '\OptIncData.mat'],'optdInc','stdInc')