%% prop ECI Two-Body: for vs. vectorized vs. for2
WC = WalkerConstellation(1,1,0,55,1000);
P = Propagator(WC);
T = 0:1:86400;
tic
[~,X] = P.prop_ECI_TB(T);
toc
tic
[~,X1] = P.prop_ECI_TB_for(T);
toc
sum(sqrt(dot(X-X1,X-X1,1)))
WC = WalkerConstellation();
P = Propagator(WC);
T = 0:1:86400;
tic
[~,X] = P.prop_ECI_TB(T);
toc
tic
[~,X1] = P.prop_ECI_TB_for(T);
toc
tic
[~,X2] = P.prop_ECI_TB_for2(T);
toc
sum(sqrt(dot(X-X1,X-X1,1)))
sum(sqrt(dot(X-X2,X-X2,1)))
%% prop ECI J2
WC = WalkerConstellation();
P = Propagator(WC);
T = 0:1:86400;
tic
[~,X] = P.prop_ECI_J2(T);
toc