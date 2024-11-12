clear
oeRange = [3000,0.2,179,360,360,360].';
minPer = 7000;
minEcc = 0.003;
primary = earth();
J2 = primary.J2;
%% Internal Test
nTest = 200000; % 10,000 take 2500s in parallel
fValVec = nan(1,nTest);
oeMat = nan(6,nTest);
% fValSum = 0;
% iterCount = 0;
% funcCount = 0;
tic
parfor iTest = 1:nTest
    rNum = rand(6,1);
    oe = oeRange.*rNum;
    oe(2) = oe(2) + minEcc; % minimum eccentricity - no errors for e>0.01
    % make sure perigee is out of earth
    rP = oe(1)+minPer;
    oe(1) = rP/(1-oe(2));
    
    % oe is mean
    [oeM,fVal] = me2oscNum(oe);
    fValVec(iTest) = fVal;
    oeMat(:,iTest) = oe;

    % Test Part
    % [oeOsc1,fVal,exitflag,output] = me2oscNum(oe);
    % fValSum = fValSum+fVal;
    % iterCount = iterCount + output.iterations;
    % funcCount = funcCount + output.funcCount;
    % if fVal > 1e-4
    %     break
    % end
end
toc
% meanfVal = fValSum/nTest
% meanIter = iterCount/nTest
% meanfunc = funcCount/nTest
probM = oeMat(:,fValVec>1e-4);
fValM = fValVec(fValVec>1e-4);
semilogy(fValM,'o')

%% ReAnalyse Problems
nProb = size(probM,2);
fVal2 = nan(1,nProb);
parfor iProb = 1:nProb
    [~,fVal] = me2oscNum(probM(:,iProb));
    fVal2(iProb) = fVal;
end
fValRatio = fVal2./fValM;
figure(1)
semilogy(1:nProb,fValM,'o',1:nProb,fVal2,'o')
figure(2)
semilogy(1:nProb,fValRatio,'o')