clear
oeRange = [3000,0.2,179,360,360,360].';
minPer = 7000;
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
    oe(2) = oe(2) + 0.001; % minimum eccentricity - no errors for e>0.01
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
%% Analyse Problems
nProbM = size(probM,2);
fValM = nan(1,nProbM);
diffM = nan(size(probM));
for iProb = 1:nProbM
    [oeOscM,fVal] = me2oscNum(probM(:,iProb));
    fValM(iProb) = fVal;
    diffM(:,iProb) = oeOscM - probM(:,iProb);
end
semilogy(fValM,'o')
