function xoverKids  = LatIntCrossover(parents,options,GenomeLength,...
    FitnessFcn,unused,thisPopulation)

%IntCon constraints
intCon = [1, 2, 3, 4 ,5];
range = options.PopInitRange;
% How many children to produce?
nKids = length(parents)/2;
% Allocate space for the kids
xoverKids = zeros(nKids,GenomeLength);
% To move through the parents twice as fast as the kids are
% being produced, a separate index for the parents is needed
iParent = 1;
% for each kid...
for iKid = 1:nKids
    % get parents
    r1 = parents(iParent);
    iParent = iParent + 1;
    r2 = parents(iParent);
    iParent = iParent + 1;
    % Children are arithmetic mean of two parents
    % ROUND will guarantee that they are integer.
    alpha = rand;
    xoverKids(iKid,:) = alpha*thisPopulation(r1,:) + ...
        (1-alpha)*thisPopulation(r2,:);
    xoverKids(iKid,intCon) = round(xoverKids(iKid,intCon));
    % Verify nPlanes is a divisor of nSats
    pList = divisors(xoverKids(iKid,1));
    diff = pList - xoverKids(iKid,2);
    [~,iClose] = min(abs(diff));
    xoverKids(iKid,2) = xoverKids(iKid,2) + diff(iClose);
    % Shift Phasing parameters to correct range
    xoverKids(iKid,[3,5]) = 1 + mod(xoverKids(iKid,[3,5]),xoverKids(iKid,2));
    nAops = xoverKids(iKid,1)/xoverKids(iKid,2);
    xoverKids(iKid,4) = 1 + mod(xoverKids(iKid,4),nAops);
    % Check inc WHY?!?
    if xoverKids(iKid,6) < range(1,6)
        xoverKids(iKid,6) = range(1,6);
    elseif xoverKids(iKid,6) > range(2,6)
        xoverKids(iKid,6) = range(2,6);
    end
    % Check that hA is above hA circ
    primary = earth();
    hACirc = CalcRgtSmaApoHeight(xoverKids(iKid,6),0,14,1)-primary.Re;
    if xoverKids(iKid,7) < hACirc
        xoverKids(iKid,7) = hACirc + 1e-3;
    end
end