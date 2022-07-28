clear
%% Test parameters
nTime = 300;
nOrb = 1;
nErr = 0;
nTest = 100000;
maxSma = 25000;

hErr = inf(nTest,1);

%% Physical Parameters
primary = earth();
mu = primary.mu;
Re = primary.Re;
J2 = primary.J2;
dScale = 1;
tScale = 1;

%% Test
tic
for iTest = 1:nTest
    r = rand(5,1);
    ecc = r(2)*0.7+0.001;
    minSma = (primary.Re+100)/(1-ecc);
    sma = minSma + r(1)*(maxSma-minSma);
    inc = r(3)*180;
    ran = r(4)*360;
    aop = r(5)*360;
    man = 0;
    f = me2ta(man,ecc);
    
    % Coordinate Switch
    radQ = sma*((1-ecc^2)/(1+ecc*cosd(f)));   % r
    aolQ = pi/180*(aop + f);                  % theta
    ranQ = pi/180*ran;                        % nu
    vraP = sqrt(mu/sma/(1-ecc^2))*ecc*sind(f);% R
    amoP = sqrt((1-ecc^2)*sma*mu);            % Theta
    amzP = amoP*cosd(inc);                    % 
    h = 0.5*vraP^2 + 0.5*amoP^2/radQ^2 - mu/radQ + ...
    0.25*mu*J2*Re^2/radQ^3*(1-3*amzP^2/amoP^2);

    % Propagate
    IC = [sma,ecc,inc,ran,aop,man];
    Sat = SingleSat(IC);
    Prop = Propagator(Sat);
    try
        [~,oe,hVec] = Prop.PropOeDeprit(nTime,nOrb);
        hErr(iTest) = max(abs(hVec-h));
    catch 
        nErr = nErr + 1;
    end
    
end
toc