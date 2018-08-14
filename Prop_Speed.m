%% prop ECI Two-Body: for vs. vectorized vs. for2
WC = WalkerConstellation(1,1,0,55,1000);
P = Propagator(WC);
T = 0:1:86400;
tic
[~,X] = P.prop_ECI_TB(T);
toc
tic
[~,X1] = P.prop_ECI_TB_for2(T);
toc
sum(sqrt(dot(X-X1,X-X1,1)));
WC = WalkerConstellation(300,10,2,55,1000);
P = Propagator(WC);
T = 0:1:86400;
tic
[~,X] = P.prop_ECI_TB(T);
toc
% tic
% [~,X1] = P.prop_ECI_TB_for(T);
% toc
tic
[~,X2] = P.prop_ECI_TB_for2(T);
toc
sum(sqrt(dot(X-X1,X-X1,1)));
sum(sqrt(dot(X-X2,X-X2,1)));
%% prop ECI J2
WC = WalkerConstellation();
P = Propagator(WC);
T = 0:1:86400;
disp('J2 GPS:')
tic
[Time,X] = P.PropEciJ2(T);
toc
WC = WalkerConstellation(200,20,2,55,1000);
P = Propagator(WC,1e-6,1e-6);
T = 0:1:60*60*2;
disp('J2 Large Con:')
tic
[Time,X] = P.PropEciJ2(T);
toc
%% prop OE Mean
% WC = SingleSat([10000,0.7,30,100,50,0]);
WC = WalkerConstellation(40,20,5,55,1000);
P = Propagator(WC);
T = 0:1:86400;
disp('Mean ode:')
tic
[~,X] = P.prop_OE_Mean(T);
toc
disp('Mean Linear:')
tic
[~,X1] = P.prop_OE_Mean_lin(T);
toc
disp(['Err = ' num2str(norm(X(:,1:6)-X1(:,1:6)))])

[~,X_ECI] = P.prop_ECI_J2(T);
M = size(X_ECI,1);
X_ECI = reshape(X_ECI.',6,P.Con.N_sats*M);
X_ECI = real(eci2oe(X_ECI(1:3,:),X_ECI(4:6,:),P.Con.mu));
X_ECI(6,:) = ta2me(X_ECI(6,:),X_ECI(2,:));
X_ECI = osc2me(X_ECI);
X_ECI = reshape(X_ECI,6*P.Con.N_sats,M).';
plot(T,wrapTo360(X1(:,5)+X1(:,6)),T,wrapTo360(X_ECI(:,5)+X_ECI(:,6)))
% legend('O1','w1','M1','O2','w2','M2')
% plot(T,X(:,4:6))
%% prop OE Osc vs prop ECI
WC = WalkerConstellation(6,3,1,55,1000);
P = Propagator(WC);
T = [0 86400];
tic
[ECI_T,X_ECI] = P.prop_ECI_J2(T);
eci_time = toc;
M = size(X_ECI,1);
X_ECI = reshape(X_ECI.',6,P.Con.N_sats*M);
X_ECI = eci2oe(X_ECI(1:3,:),X_ECI(4:6,:),P.Con.mu);
X_ECI(6,:) = ta2me(X_ECI(6,:),X_ECI(2,:));
X_ECI = reshape(X_ECI,6*P.Con.N_sats,M).';
tic
[oe_T,X_OE] = P.prop_OE_Osc(T);
oe_time = toc;
diff = sum(sqrt(dot(X_ECI-X_OE,X_ECI-X_OE,1)));
disp([newline 'ECI Time: ' num2str(eci_time) newline ' OE Time: ' num2str(oe_time) newline ' Error: ' num2str(diff)])
plot(T,X_ECI(:,3:6),T,X_OE(:,3:6))