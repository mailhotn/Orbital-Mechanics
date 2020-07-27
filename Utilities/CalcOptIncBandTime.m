function inc = CalcOptIncBandTime(latEm, delLat, elevMin)
lMin = latEm - delLat;
lMax = latEm + delLat;
inc = fminbnd(@(inc) mincost(inc,lMax,lMin, elevMin), latEm, 90);
end

function cost = mincost(inc,lMax,lMin,elevMin)
sma = CalcRgtSma(0,inc,14,1);
primary = earth();
H = sma-primary.Re;
rho = asind(primary.Re/(primary.Re+H));
eta = asind(cosd(elevMin)*sind(rho));
lam = 90 - elevMin - eta;

cost = abs(inc - (lMin+lam));
% cost = -(asind(sind(lMax+lam)/sind(inc))-asind(sind(lMin-lam)/sind(inc)));
end