function [c,ceq] = LatticeGaCollCon(x,Arch,latEm,OptParams)
%LatticeGaCollCon Checks for collisions in lattice for GA
ceq = [];
% Handle Input
hA = OptParams.hAList(x(1));
inc = latEm + x(2);
Phase = struct();
Phase.nC1 = x(3);
Phase.nC2 = x(4);
Phase.nC3 = x(5);

% Orbital Parameters
[sma,ecc] = CalcRgtSmaApoHeight(inc, hA, OptParams.nRepeats, OptParams.nDays);
Orbit = struct();
Orbit.sma = sma;
Orbit.ecc = ecc;
Orbit.inc = inc;
Orbit.hA = hA;

% Initial Conditions
InitCon = InitConElliptical(ecc,inc,sma,latEm,0);

% Creat Constellation
LC = LatticeConstellation(Arch,Phase,Orbit,InitCon);

% Calculate Distance
minDist = CalcMinDist(LC);

c = 1 - minDist;

end

