function Cz = StumC(z)
% StumC calculates the C Stumpff function of input z
if z > 0
    Cz = (1-cos(sqrt(z)))/z;
elseif z < 0
    Cz = (cosh(sqrt(-z))-1)/(-z);
else
    Cz = 1/2;
end
end