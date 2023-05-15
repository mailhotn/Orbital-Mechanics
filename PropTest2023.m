clear
primary = earth();
mu = primary.mu;
J2 = primary.J2;
Re = primary.Re;

%% Define Large Constellation
Arch.nPlanes = 6;
Arch.nAops = 4;
Arch.nSatsPerAop = 1;

Phase.nC1 = 5;
Phase.nC2 = 3;
Phase.nC3 = 6;

Orbit.sma = 7000+rand*10000;
Orbit.ecc = 0.001 + rand*0.5;
Orbit.inc = 90*rand;

Con = LatticeConstellation(Arch,Phase,Orbit);
Prop = Propagator(Con);

%% Propagate
t = 0:10:86400;
tic
[~,eci0] = Prop.PropEciJ2(t);
toc
eci0 = reshape(eci0.',6,numel(eci0)/6);
oe0 = eci2oe(eci0);
tic
[~,oe1] = Prop.PropOeOsc(t);
toc
oe1 = reshape(oe1.',6,numel(oe1)/6);
oe1(6,:) = wrapTo360(me2ta(oe1(6,:),oe1(2,:)));
eci1 = oe2eci(oe1);
tic
[~,oe2] = Prop.PropOeOsc3(t);
toc
oe2 = reshape(oe2.',6,numel(oe2)/6);
oe2(6,:) = wrapTo360(me2ta(oe2(6,:),oe2(2,:)));
eci2 = oe2eci(oe2);
tic
[~,oe3, pns] = Prop.PropOePns(t);
toc
pns = reshape(pns.',6,numel(pns)/6); 
oe3 = reshape(oe3.',6,numel(oe3)/6);
oe3 = me2ta(oe3);
%% Error Analysis
errMax1 = max(abs(oe1-oe0),[],2) %#ok<*NOPTS>
errMen1 = mean(abs(oe1-oe0),2)

errMax2 = max(abs(oe2-oe0),[],2)
errMen2 = mean(abs(oe2-oe0),2)

errMax3 = max(abs(oe3-oe0),[],2)
errMen3 = mean(abs(oe3-oe0),2)

errMax12 = max(abs(oe1-oe2),[],2)
errMen12 = mean(abs(oe1-oe2),2)

errMax13 = max(abs(oe1-oe3),[],2)
errMen13 = mean(abs(oe1-oe3),2)

%% Error by Satellite

%% Hamiltonian Analysis
H0 = -mu./(2*oe0(1,:))+mu*J2*Re^2./(2*oe0(1,:).^3.*(1-oe0(2,:).^2).^3).*...
    (1+oe0(2,:).*cosd(oe0(6,:))).^3.*...
    (3*sind(oe0(3,:)).^2.*sind(oe0(6,:)+oe0(5,:)).^2-1);
H0 = reshape(H0,Con.nSats,length(t));

N0 = cosd(oe0(3,:)).*sqrt(mu*oe0(1,:).*(1-oe0(2,:).^2));
N0 = reshape(N0,Con.nSats,length(t));

H1 = -mu./(2*oe1(1,:))+mu*J2*Re^2./(2*oe1(1,:).^3.*(1-oe1(2,:).^2).^3).*...
    (1+oe1(2,:).*cosd(oe1(6,:))).^3.*...
    (3*sind(oe1(3,:)).^2.*sind(oe1(6,:)+oe1(5,:)).^2-1);
H1 = reshape(H1,Con.nSats,length(t));

N1 = cosd(oe1(3,:)).*sqrt(mu*oe1(1,:).*(1-oe1(2,:).^2));
N1 = reshape(N1,Con.nSats,length(t));

H2 = -mu./(2*oe2(1,:))+mu*J2*Re^2./(2*oe2(1,:).^3.*(1-oe2(2,:).^2).^3).*...
    (1+oe2(2,:).*cosd(oe2(6,:))).^3.*...
    (3*sind(oe2(3,:)).^2.*sind(oe2(6,:)+oe2(5,:)).^2-1);
H2 = reshape(H2,Con.nSats,length(t));

N2 = cosd(oe2(3,:)).*sqrt(mu*oe2(1,:).*(1-oe2(2,:).^2));
N2 = reshape(N2,Con.nSats,length(t));

H3 = 0.5*(pns(4,:).^2 +pns(5,:).^2./pns(1,:).^2) -mu./pns(1,:) + ...
    0.5*mu*J2*Re^2./pns(1,:).^3.*...
    (3*(pns(6,:).^2./pns(5,:)-pns(6,:).^4./(4*pns(5,:).^2)).*...
    sin(pns(2,:)+pns(3,:)./pns(6,:)).^2-1);
H3 = reshape(H3,Con.nSats,length(t));

N3 = pns(5,:)-pns(6,:).^2/2;
N3 = reshape(N3,Con.nSats,length(t));