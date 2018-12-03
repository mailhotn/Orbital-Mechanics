function [ ecc ] = CalcFrozenEcc( sma, inc, aop )
%CalcFrozenEcc calculates e for a frozen orbit

%   ~~~~~~~~~~~~~  INPUTS  ~~~~~~~~~~~~~~
%   sma   - semimajor axis [km]
%   inc   - orbital inclination [deg]
%   aop   - argument of perigee [deg]

% Handle Input
inc = inc*pi/180;
aop = aop*pi/180;
% Planet
primary = earth();
J2 = primary.J2;
J3 = primary.J3;
Re = primary.Re;
% N-R intialization
tol = 1e-10;
ratio = inf;
iter = 0;
maxIter = 20;

x = J2;
while abs(ratio) > tol && iter < maxIter
    f = J2 + J3/2*Re*sin(aop)/sma*...
        1/(x - x^3)*...
        (sin(inc) - x*cos(inc)^2/sin(inc));
    df = J3/2*Re*sin(aop)/sma*(...
        (3*x^2 - 1)/(x - x^3)^2*sin(inc) -...
        2*x/(1 - x^2)^2 *cos(inc)^2/sin(inc));
    ratio = f/df;
    x = x - ratio;
    iter = iter + 1;
end

if imag(x) == 0 && real(x) > 0
    ecc = x;
else
    error(['No viable solution found! inc = ' num2str(inc*180/pi) 'sma = ' num2str(sma)])
end