function [ E ] = M2E( M , e, tol)
%M2E Calculates the Eccentric anomaly (in radians) corresponding 
% to the given Mean anomaly (in radians)
maxIter = 20;
if nargin < 3
    tol = 1e-14;
end
E = (M < pi).*(M+e) + (M >= pi).*(M-e);
dE = inf(1,length(M));
iter = 0;
while (norm(dE,inf) > tol) && (iter < maxIter)
    dE = -(E-e.*sin(E)-M)./(1-e.*cos(E));
    E  = E + dE;
    iter = iter + 1;
end

E = wrapTo2Pi(E);
if iter >= maxIter
    warning(['M2E maximum iterations reached, result may not have converged'...
             newline 'Final increment dE = ' num2str(max(abs(dE)))])
end
end