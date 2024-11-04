clear
oeRange = [3000,0.2,180,360,360,360].';
minPer = 7000;
%% Internal Test
nTest = 10000;
errMatM = nan(6,nTest);
errMatO = nan(6,nTest);
% fValSum = 0;
% iterCount = 0;
% funcCount = 0;
tic
for iTest = 1:nTest
    rNum = rand(6,1);
    oe = oeRange.*rNum;
    % make sure perigee is out of earth
    rP = oe(1);
    oe(1) = rP/(1-oe(2));
    
    
    % oe is mean
    oeM = osc2meNum(me2oscNum(oe));
    errMatM(:,iTest) = oeM-oe;

    % oe is osc
    oeOsc = me2oscNum(osc2meNum(oe));
    errMatO(:,iTest) = oeOsc-oe;
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