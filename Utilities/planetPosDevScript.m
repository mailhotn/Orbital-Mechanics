t0 = juliandate(2000,1,1);

t = (t0:0.25:(t0+365)).';
ts = (t-t0);

[rM, vM] = planetEphemeris(t,'Earth','Moon');
xM = rM(:,1);
yM = rM(:,2);
zM = rM(:,3);
dxM = vM(:,1);
dyM = vM(:,2);
dzM = vM(:,3);
oeM = eci2oe3b(rM.',vM.',Earth,Moon);
