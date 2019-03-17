function [ sma, ecc, fVal] = CalcRgtSmaApoHeight( inc, hA, j, k )
%CalcRgtSmaApoHeight calculates e & a for an RGT Orbit with given apogee
%height
% Setting hA to 0 will result in a circular orbit
%   Uses multivariate Newton - Raphson method

%   ~~~~~~~~~~~~~  INPUTS  ~~~~~~~~~~~~~~
%   inc - orbital inclination [deg]
%   hA  - desired apogee height [km]
%   j   - number of repeats
%   k   - number of days to repeat

if hA == 0
    ecc = 0;
    sma = CalcRgtSma(ecc,inc,j,k);
    fVal = 0;
    return
end
% Planet
primary = earth();
J2 = primary.J2;
Re = primary.Re;
mu = primary.mu;
Tday = primary.Tday;
we = primary.we;

% Handle Input
inc = inc*pi/180;
rA = hA + primary.Re;

% N-R intialization
tol = 1e-12;
delX = inf(2,1);
fVal = ones(2,1);
iter = 0;
maxIter = 20;

initSma = (((k*Tday/j)/2/pi)^2*mu)^(1/3);
initEcc = rA/initSma - 1;
x = [initEcc;
    initSma];

while norm(delX,inf) > tol && iter < maxIter
    J2Effects = 3/4*J2*Re^2/(1 - x(1)^2)^2 * ...
        (2*j/k*cos(inc) - sqrt(1 - x(1)^2)*(3*cos(inc)^2 - 1) - (5*cos(inc)^2 - 1));
    
    f1 = rA - x(2)*(1 + x(1));
    
    f2 = x(2)^2 - j*we*(pi/180)/(k*sqrt(mu))*x(2)^(7/2) - J2Effects;
    
    df1de = -x(2);
    
    df1da = -(1 + x(1));
    
    df2de = 3/4*J2*Re^2*...
        x(1)/(1 - x(1)^2)^3*...
        (-4*(5*cos(inc)^2 - 1) + 8*j/k*cos(inc) - sqrt(1-x(1)^2)*...
        3*(3*cos(inc)^2 - 1));
    
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