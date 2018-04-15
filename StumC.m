function Cz = StumC(z,a)
% StumC calculates the C Stumpff function of input z
% In order for this to work for symbolic inputs when solving the universal
% Kepler equation, the sign is taken from alpha = 1/a, the inverse of the
% semimajor axis, determining the type of orbit
if nargin == 1
    a = z;
end
if a > 0
    Cz = (1-cos(sqrt(z)))/z;
elseif a < 0
    Cz = (cosh(sqrt(-z))-1)/(-z);
else
    Cz = 1/2;
end
end