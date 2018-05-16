function [c, ceq] = Walker_GA_nonlcon(x, T)
%Walker_GA_nonlcon Enforces mod(T,P)==0
%   T/P must be 0 for a valid Walker Constellation.
%   c <=0 is used, since ceq cannot be used if the genome x includes
%   integer values. This is valid since mod(T,P) >= 0.
ceq = [];
P = x(1);
F = x(2);
c(1) = mod(T,P);
c(2) = (F-P+1);
end

