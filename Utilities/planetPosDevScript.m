primary = Earth;
third = Sun;
t0 = juliandate(2000,1,1);

t = (t0:0.25:(t0+365)).';
ts = (t-t0);

[rM, vM] = planetEphemeris(t,primary.name,third.name);
xM = rM(:,1);
yM = rM(:,2);
zM = rM(:,3);
dxM = vM(:,1);
dyM = vM(:,2);
dzM = vM(:,3);
[rFit, vFit] = third.PosJ2000(ts*86400);
[oeM,EM] = eci2oe3b(rM,vM,primary,third);
[oeFit,EFit] = eci2oe3b(rFit,vFit,primary,third);

%% 
figure(1)
tiledlayout(3,2,"Padding","compact","Tilespacing","tight")
nexttile
plot(t,oeM(1,:),t,oeFit(1,:))
nexttile
plot(t,oeM(2,:),t,oeFit(2,:))
nexttile
plot(t,oeM(3,:),t,oeFit(3,:))
nexttile
plot(t,unwrap(oeM(4,:)),t,unwrap(oeFit(4,:)))
nexttile
plot(t,oeM(5,:),t,oeFit(5,:))
nexttile
plot(t,unwrap(oeM(6,:)+oeM(5,:)),t,unwrap(oeFit(6,:)+oeFit(5,:)))