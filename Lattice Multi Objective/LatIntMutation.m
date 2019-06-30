function mutants = LatIntMutation(parents, options, GenomeLength, ...
    ~, state, ~, thisPopulation)

% Function that creates the mutated children using the Gaussian
% distribution. It does not satisfy linear constraints!

%IntCon constraints
IntCon = [1, 2, 3, 4 ,5];

rate = 0.05;
shrink = 0.5; 
scale = 1;
scale = scale - shrink * scale * state.Generation/options.Generations;

mutationPop =  length(parents);
mutants =  zeros(mutationPop, GenomeLength);

indMut = rand(mutationPop,5) < rate;
for iMut = 1:mutationPop
    parent = thisPopulation(parents(iMut),:);
    randNum = randn(1,length(parent));
    % Mutate nSats - Creep
    if indMut(iMut,1)
    mutants(iMut,1) = round(parent(1) + scale*5*randNum(1));
    else 
    % Mutate nPlanes
    pList = divisors(mutants(iMut,1));
    mutants(iMut,2) = round(parent(2) + scale*(pList(end) - 1)*randn);
    % Verify nPlanes is a divisor of nSats
    diff = pList - mutants(iMut,2);
    [~,iClose] = min(abs(diff));
    mutants(iMut,2) = mutants(iMut,2) + diff(iClose);
    % mutate phasing
    mutants(iMut,[3,5]) = round(parent([3,5]) + scale*mutants(iMut,2)*randNum([3,5]));
    % Shift Phasing parameters to correct range
    mutants(iMut,[3,5]) = 1 + mod(mutants(iMut,[3,5]),mutants(iMut,2));
    nAops = mutants(iMut,1)/mutants(iMut,2);
    mutants(iMut,4) = round(parent(4) + scale*nAops*randNum(4));
    mutants(iMut,4) = 1 + mod(mutants(iMut,4),nAops);
    % Mutate inc
    mutants(iMut,6) = parent(6) + scale*(upper(6)-lower(6))*randNum(6);
    % Check that hA is above hA circ
    primary = earth();
    hACirc = CalcRgtSmaApoHeight(mutants(iMut,6),0,14,1)-primary.Re;
    mutants(iMut,7) = parent(7) + scale*(1000 - hACirc)*randNum(7);
    if mutants(iMut,7) < hACirc
        mutants(iMut,7) = hACirc;
    elseif mutants(iMut,7) > 1000
        mutants(iMut,7) = 1000;
    end
end
    
