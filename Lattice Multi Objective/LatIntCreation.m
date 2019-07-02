function Population = LatIntCreation(GenomeLength, ~, options)
% LatIntCreation Function that creates an initial population satisfying bounds and
% integer constraints for Lattice Flower Constellation multi-objective GA

totalPopulation = sum(options.PopulationSize);


range = options.PopInitRange;
lower = range(1,:);
span =  range(2,:) - lower;

Population = zeros(totalPopulation,GenomeLength);

% Set nSats
Population(:,1) = repmat(lower(1),totalPopulation,1) + ...
    randi(span(1) + 1,totalPopulation,1) - 1;
for iPop = 1:totalPopulation
    % Set nPlanes
    pList = divisors(Population(iPop,1));
    Population(iPop,2) = pList(randi(numel(pList)));
    % Set Phasing
    Population(iPop,[3,5]) = randi(Population(iPop,2),1,2);
    nAops = Population(iPop,1)/Population(iPop,2);
    Population(iPop,4) = randi(nAops);
    % Set inc
    Population(iPop,6) = lower(6) + span(6)*rand;
    % Set hA
    primary = earth();
    hACirc = CalcRgtSmaApoHeight(Population(iPop,6),0,14,1)-primary.Re;
    Population(iPop,7) = hACirc + rand*(1000-hACirc - 1e-3) + 1e-3;
end