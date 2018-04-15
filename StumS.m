function Sz = StumS(z,a)
% StumS calculates the S Stumpff function of input z
% In order for this to work for symbolic inputs when solving the universal
% Kepler equation, the sign is taken from alpha = 1/a, the inverse of the
% semimajor axis, determining the type of orbit
if nargin == 1
    a = z;
end
if a > 0
    Sz = (sqrt(z)-sin(sqrt(z)))/sqrt(z)^3;
elseif a < 0
    Sz = (sinh(sqrt(-z))-sqrt(-z))/sqrt(-z)^3;
else
    Sz = 1/6;
end
end