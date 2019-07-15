datafolder = 'C:\Users\User\Dropbox\Lattice Optimization Data\Multi Objective';
intTarget = 10;
p90Target = 5;
covTarget = 99;
maxTarget = 10;
latList = 60;
%% Loop
for iLat = 1:numel(latList)
    load([datafolder '\LatticeGaSol_Lat_' num2str(latList(iLat)) '.mat'])
    achieveInt = mean(GaSol.intPdop,1) < intTarget;
    achievep90 = mean(GaSol.p90,1) < p90Target;
    achieveCov = mean(GaSol.coverage,1) > covTarget;
    achieveMax = mean(GaSol.maxPdop,1) < maxTarget;
    achieveAll = achieveInt & achievep90 & achieveCov & achieveMax;
    conList = find(achieveAll);
    
    parSats = GaSol.Cons{conList(1)}.nSats;
    parPlanes = GaSol.Cons{conList(1)}.nPlanes;
    for iCon = 2:numel(conList) % ascending nSats
        if GaSol.Cons{conList(iCon)}.nPlanes < min(parPlanes)
            parSats = [parSats GaSol.Cons{conList(iCon)}.nSats];
            parPlanes = [parPlanes GaSol.Cons{conList(iCon)}.nPlanes];
        end
    end
    figure(10*iLat + 1)
    plot(parSats,parPlanes,'--o')
end

