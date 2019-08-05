function mutants = LatIntMutation(parents, options, GenomeLength, ...
    FitnessFcn, state, score, thisPopulation)

%IntCon constraints
IntCon = [1, 2, 3, 4 ,5];
range = options.PopInitRange;
mutRate = 0.8;
shrink = 0.2;
scale = 1;
scale = scale - shrink * scale * state.Generation/options.Generations;

nSatsStd = 5; % Satellites
incStd = 5; % deg
hAStd = 50; % km

nMutants =  length(parents);
mutants =  zeros(nMutants, GenomeLength);
% array indicating which genes should mutate in which parents
genes2Mut = rand(nMutants,5) < mutRate;
randNum = randn(nMutants,3);
for iMut = 1:nMutants
    parent = thisPopulation(parents(iMut),:);
    
    %% Mutate nSats - Creep sometimes
    if genes2Mut(iMut,1) % mutate
        mutants(iMut,1) = round(parent(1) + scale*nSatsStd*randNum(iMut,1));
        if mutants(iMut,1) < range(1,1)
            mutants(iMut,1) = range(1,1);
        elseif mutants(iMut,1) > range(2,1)
            mutants(iMut,1) = range(2,1);
        end
    else % inherit
        mutants(iMut,1) = parent(1);
    end
    
    %% Mutate nPlanes - Random Reset
    pList = divisors(mutants(iMut,1));
    if genes2Mut(iMut,2) % mutate
        mutants(iMut,2) = pList(randi(numel(pList)));
    else % inherit
        mutants(iMut,2) = parent(2);
    end
    % Verify nPlanes is a divisor of nSats
    % This fixes nPlanes if nSats mutated and nPlanes is no longer a divisor
    diff = pList - mutants(iMut,2);
    [~,iClose] = min(abs(diff));
    mutants(iMut,2) = mutants(iMut,2) + diff(iClose);
    
    %% Mutate phasing - Random Reset
    % Mutate nC1, nC3
    mutants(iMut,[3,5]) = [(genes2Mut(iMut,3))*randi(mutants(iMut,2)) + ...% mutate
                           ~(genes2Mut(iMut,3))*parent(3),...% inherit
                           (genes2Mut(iMut,5))*randi(mutants(iMut,2)) + ...% mutate
                           ~(genes2Mut(iMut,5))*parent(5)];% inherit
    % Shift Phasing parameters to correct range
    % Also fixes phasing if nPlanes changed
    mutants(iMut,[3,5]) = 1 + mod(mutants(iMut,[3,5]),mutants(iMut,2));
    
    % Mutate nC2
    nAops = mutants(iMut,1)/mutants(iMut,2);
    mutants(iMut,4) = (genes2Mut(iMut,4))*randi(nAops) + ...% mutate
                           ~(genes2Mut(iMut,4))*parent(4);% inherit
    % Shift to range
    mutants(iMut,4) = 1 + mod(mutants(iMut,4),nAops);
    
    %% Mutate inc - Creep always
    mutants(iMut,6) = parent(6) + scale*incStd*randNum(iMut,2);
    % Verify that new inc is in bounds
    if mutants(iMut,6) < range(1,6)
        mutants(iMut,6) = range(1,6);
    elseif mutants(iMut,6) > range(2,6)
        mutants(iMut,6) = range(2,6);
    end
    %% Mutate hA - Creep always
    mutants(iMut,7) = parent(7) + scale*hAStd*randNum(iMut,3);
    % Check that hA is above hA circ
    primary = earth();
    hACirc = CalcRgtSmaApoHeight(mutants(iMut,6),0,14,1)-primary.Re;
    if mutants(iMut,7) < hACirc
        mutants(iMut,7) = hACirc + 1e-3; % Circular orbit + 1 meter to avoid numerical issues
    elseif mutants(iMut,7) > 1000
        mutants(iMut,7) = 1000;
    end
end

