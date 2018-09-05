function Sz = StumS(z)
% StumS calculates the S Stumpff function of input z
if z > 0
    Sz = (sqrt(z)-sin(sqrt(z)))/sqrt(z)^3;
elseif z < 0
    Sz = (sinh(sqrt(-z))-sqrt(-z))/sqrt(-z)^3;
else
    Sz = 1/6;
end
end