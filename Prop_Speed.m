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
disp('J2 GPS:')
tic
[~,~] = P.prop_ECI_J2(T);
toc
WC = WalkerConstellation(1000,20,5,55,1000);
P = Propagator(WC);
T = 0:1:86400;
disp('J2 Large Con:')
tic
[~,~] = P.prop_ECI_J2(T);
toc
%% prop OE Mean
WC = WalkerConstellation(100,20,5,55,500);
P = Propagator(WC);
T = 0:60:86400*30;
disp('Mean ode:')
tic
[~,X] = P.prop_OE_Mean(T);
toc
disp('Mean Linear:')
tic
[~,X1] = P.prop_OE_Mean_lin(T);
toc
disp(['Err = ' num2str(norm(X(:,1:6)-X1(:,1:6)))])
plot(T,X(:,4:6))