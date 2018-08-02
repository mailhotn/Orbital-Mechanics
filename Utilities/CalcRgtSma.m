function [ x ] = CalcRgtSma(e, inc, j, k, primary , tol)
%CalcRgtSma calculates the semimajor axis for an RGT orbit
%   ~~~~~~~~~~~~~  INPUTS  ~~~~~~~~~~~~~~
%   e     - eccentricity
%   inc   - orbital inclination [deg]
%   j     - number of repeats
%   k     - number of days to repeat

if nargin < 5 % default is Earth
    primary = earth();
end
if nargin < 6 % default is Earth
    tol = 1e-12;
end
ratio = inf;
iter = 0;
maxIter = 20;
% N-R Initializations
J2Effects = 3/4*primary.J2*primary.Re^2/(1-e^2)^2 * ...
    (2*j/k*cosd(inc) - sqrt(1-e^2)*(3*cosd(inc)^2 - 1) - (5*cosd(inc)^2 - 1));

a = (((k*primary.Tday/j)/2/pi)^2*primary.mu)^(1/3);
while abs(ratio) > tol && iter < maxIter % Newton - Raphson
    f  = a^2 - j*primary.we*(pi/180)/(k*sqrt(primary.mu))*a^(7/2) - J2Effects;
    df = 2*a - 7*j*primary.we*(pi/180)/(2*k*sqrt(primary.mu))*a^(5/2);
    ratio = f/df;
    a = a - ratio;
    iter = iter + 1;
end
% Output Value
if imag(a) == 0
    x = a;
else
    error(['No viable solution found! inc = ' num2str(inc) 'j = ' num2str(j)])
end
end

