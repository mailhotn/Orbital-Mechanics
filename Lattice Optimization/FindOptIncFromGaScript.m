% Useless before restoring inc data
clear
datafolder =     'C:\Users\User\Dropbox\Lattice Optimization Data\GA Standard\Previous Runs\Version 2 - all lats';
load([datafolder '\OptParams']);
%%
for iLat = 1:numel(latList)
    incList = nan(1,maxSats-minSats+1);
    haList = nan(1,maxSats-minSats+1);
    
    for nSats = minSats:maxSats
        load([datafolder '\LatticeGaSol_Lat_' num2str(latList(iLat)) '_nSats_' ...
            num2str(nSats) '.mat']);
        incList(nSats-minSats+1) = GaSol.Orbit.inc;
        haList(nSats-minSats+1) = GaSol.Orbit.hA;
    end
    figure(iLat*10+1)
    plot(minSats:maxSats,incList-latList(iLat),'o')
    
    figure(iLat*10+2)
    plot(minSats:maxSats,haList,'o')
    mean(incList-latList(iLat))
    std(incList-latList(iLat))
end
    