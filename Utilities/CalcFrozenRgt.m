function [ sma, ecc, fVal] = CalcFrozenRgt( inc, aop, j, k )
%CalcFrozenRgt calculates e & a for a frozen RGT Orbit
%   Uses multivariate Newton - Raphson method

%   ~~~~~~~~~~~~~  INPUTS  ~~~~~~~~~~~~~~
%   inc - orbital inclination [deg]
%   aop - argument of perigee [deg]
%   j   - number of repeats
%   k   - number of days to repeat

% Handle Input
inc = inc*pi/180;
aop = aop*pi/180;
% Planet
primary = earth();
J2 = primary.J2;
J3 = primary.J3;
Re = primary.Re;
mu = primary.mu;
Tday = primary.Tday;
we = primary.we;

% N-R intialization
tol = 1e-12;
delX = inf(2,1);
fVal = ones(2,1);
iter = 0;
maxIter = 20;

initSma = (((k*Tday/j)/2/pi)^2*mu)^(1/3);
initEcc = -J3/2/J2*Re/initSma*sin(inc);
x = [initEcc;
    initSma];

while norm(delX) > tol && iter < maxIter
    J2Effects = 3/4*J2*Re^2/(1 - x(1)^2)^2 * ...
        (2*j/k*cos(inc) - sqrt(1 - x(1)^2)*(3*cos(inc)^2 - 1) - (5*cos(inc)^2 - 1));
    
    f1 = J2 + J3/2*Re*sin(aop)/x(2)*...
        1/(x(1) - x(1)^3)*...
        (sin(inc) - x(1)*cos(inc)^2/sin(inc));
    
    f2 = x(2)^2 - j*we*(pi/180)/(k*sqrt(mu))*x(2)^(7/2) - J2Effects;
    
    df1de = J3/2*Re*sin(aop)/x(2)*(...
        (3*x(1)^2 - 1)/(x(1) - x(1)^3)^2*sin(inc) -...
        2*x(1)/(1 - x(1)^2)^2 *cos(inc)^2/sin(inc));
    
    df1da = -J3/2*Re*sin(aop)/x(2)^2*...
        1/(x(1) - x(1)^3)*...
        (sin(inc) - x(1)*cos(inc)^2/sin(inc));
    
    df2de = -3/4*J2*Re^2*...
        x(1)/(1 - x(1)^2)^2*...
        (8*(5*cos(inc)^2 - 1) - 16*j/k*cos(inc) +9/2*(3*cos(inc)^2 - 1));
    
    df2da = 2*x(2) - 7*j*we*(pi/180)/(2*k*sqrt(mu))*x(2)^(5/2);
    
    J = [df1de df1da;
        df2de df2da];
    fVal = [f1; f2];
    delX = -J\fVal;
    
    x = x + delX;
    iter = iter + 1;
end
if x(1) > 0 && isreal(x)
    ecc = x(1);
    sma = x(2);
else
    error('No viable solution found!')
end
end