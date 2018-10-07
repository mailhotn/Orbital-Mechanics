datafolder = ['C:\Users\User\Dropbox\Lattice Optimization Data']; %#ok<NBRAK>
load([datafolder '\OptParams.mat']);

%%
ecc = 0.05;
nSats = 51;
nPlanes = 51;
latGs = 50;
load([datafolder '\LatticeExSol_Lat_' num2str(latGs) ...
    '_nSats_' num2str(nSats) '_ecc_' num2str(ecc) '.mat']);
Arch.nSats = nSats;
Arch.nPlanes = nPlanes;
iPlane = ExSol.archMat(1,:) == nPlanes;
Arch.nAops = ExSol.archMat(2,iPlane);
Arch.nSatsPerAop = ExSol.archMat(3,iPlane);

Phase.nC1 = ExSol.phaseMat(1,iPlane);
Phase.nC2 = ExSol.phaseMat(2,iPlane);
Phase.nC3 = ExSol.phaseMat(3,iPlane);

LC = LatticeConstellation(Arch,Phase,ExSol.Orbit,ExSol.InitCon);
pOrbit = LC.sma*(1-LC.ecc^2);
nOrbit = sqrt(LC.mu/LC.sma^3);
aopRate = 3/4*LC.J2*(LC.Re/pOrbit)^2*nOrbit*(5*cosd(LC.inc)^2 - 1);
aopTime = 2*pi/abs(aopRate)/86164;
nRelTraj = Arch.nAops*(nPlanes/gcd(nPlanes,Phase.nC3))
%%
timeVec = [];
oneDay = 0:10:86164;
dayStep = 5;
dayVec = 1:dayStep:aopTime;
intPdop = nan(1,length(dayVec));
coverage = nan(1,length(dayVec));
for iDay = 1:dayStep:dayVec(end)
    timeVec = [timeVec,oneDay + iDay*86164];
end

Prop = Propagator(LC,PropParams.relTol,PropParams.absTol);
tic
[Time, meanState] = Prop.PropOeMeanFast(timeVec);
toc
aop0 = meanState(:,5).';
tic
meanState = reshape(meanState.',6,length(Time)*nSats);
oscState = me2osc(meanState);
oscState(6,:) = me2ta(oscState(6,:),oscState(2,:));

[R,V] = oe2eci(oscState);
eciState = reshape([R;V],6*nSats,length(Time));
toc
tic
pdop = TdoaPdopVec(eciState.',Time,latGs,0,0,PropParams.elevMin);
toc
for iDay = 1:(length(dayVec))
    pdopDay = pdop((length(oneDay)*(iDay-1)+1):length(oneDay)*iDay);
    coverage(iDay) = sum(~isnan(pdopDay))/length(pdopDay)*100;
    pdopDay(isnan(pdopDay)) = 1000;
    intPdop(iDay) = trapz(oneDay,pdopDay)./(oneDay(end)-oneDay(1));
%     if rem(iDay-1,50) == 0
%         PlotGroundTrack(eciState(1:6,(length(oneDay)*(iDay-1)+1):length(oneDay)*iDay).',...
%             oneDay,0);
%     end
end
% plot(dayVec(1:ceil(length(dayVec)/2)),meanPdop(1:ceil(length(dayVec)/2))...
%     ,dayVec(ceil(length(dayVec)/2):end)-dayVec(ceil(length(dayVec)/2))...
%     ,meanPdop(ceil(length(dayVec)/2):end)...
%     ,'Linewidth',1.5)
figure(2)
semilogy(dayVec,intPdop,'.','Linewidth',1.5)
grid minor
xlabel('Time [Days]')
ylabel('$\frac{1}{T_d}\int^{T_d}_{0}{PDOPdt}$','interpreter','latex',...
    'fontsize',12)
figure(3)
plot(dayVec,coverage,'.','Linewidth',1.5)
grid minor
xlabel('Time [Days]')
ylabel('Coverage %')
ylim([0,100])