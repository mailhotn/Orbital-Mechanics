clear
icM = [7000,0.01,30,10,10,10].';
tOe = [7100,0.01,30].';
icOsc = me2osc(icM);
primary = Earth;
dScale = primary.Re;
tScale = sqrt(primary.Re^3/primary.mu);
fScale = dScale/tScale^2;

thrustMag = 0.00001;
tol = 1e-5;

Sat = SingleSat(icOsc,primary);
Control = ZhangInverse(tOe,thrustMag,tol);
Prop = Propagator(Sat,Control);

t = 0:10:12000;

[~,X] = Prop.PropConOeOsc(t);

