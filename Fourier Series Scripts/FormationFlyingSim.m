clear
oeLead = [7000, 0.01, 5, 10, 10, 10].'; % init osc oe
% oeFollow = f(oeLead)

leadSat = SingleSat(oeLead);
leadProp = Propagator(leadSat);

t = 0:10:86400;

%% Prop
[~,oeC] = leadProp.PropOeOsc3(t);
oeC = oeC.';
[~,oeM] = leadProp.PropOeMeanShort(t);
oeB = me2oscSP(oeM.');
[~,oeF] = leadProp.PropOeFourier2Ord(t,5);
oeF = oeF.';
errB = abs(oeC-oeB);
errF = abs(oeC-oeF);

