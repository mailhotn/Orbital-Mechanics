function [c, ceq] = Walker_GA_nonlcon(x, T)
%Walker_GA_nonlcon Enforces mod(T,P)==0
%   rem(T/P) must be 0 for a valid Walker Constellation.
%   c <=0 is used, since ceq cannot be used if the genome x includes
%   integer values. This is valid since mod(T,P) >= 0.
ceq = [];
P = abs(x(1));
F = abs(x(2));
c(1) = mod(T,P)*100; % Added gain to help with failure when mod == 1. 
                     % Unclear if this works.
c(2) = (F-P+1);
end

