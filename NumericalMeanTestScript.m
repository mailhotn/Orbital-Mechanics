clear
oeRange = [3000,0.2,180,359,359,359].';
minPer = 7000;
%% Internal Test
nTest = 100000;
errMatM = nan(6,nTest);
errMatO = nan(6,nTest);
oeMat = nan(6,nTest);
% fValSum = 0;
% iterCount = 0;
% funcCount = 0;
tic
parfor iTest = 1:nTest
    rNum = rand(6,1);
    oe = oeRange.*rNum;
    oe(2) = oe(2) + 0.01;
    % make sure perigee is out of earth
    rP = oe(1)+minPer;
    oe(1) = rP/(1-oe(2));
    
    % oe is mean
    oeM = osc2meNum(me2oscNum(oe));
    errMatM(:,iTest) = oeM-oe;

    % oe is osc
    oeOsc = me2oscNum(osc2meNum(oe));
    errMatO(:,iTest) = oeOsc-oe;
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
