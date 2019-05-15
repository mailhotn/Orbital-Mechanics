function [ sma, ecc, inc, fVal ] = RgtSunSynElements( hA, j, k )
%RgtSunSynElements calculates i, e & a for an RGT, sun-synchronous Orbit 
%with given apogee height
% Setting hA to 0 will result in a circular orbit
%   Uses multivariate Newton - Raphson method

%   ~~~~~~~~~~~~~  INPUTS  ~~~~~~~~~~~~~~
%   hA  - desired apogee height [km]
%   j   - number of repeats
%   k   - number of days to repeat

if hA == 0
    [sma, ecc, inc, fVal] = RgtSunSynCirc(j, k);
    return
end
% Planet
primary = earth();
J2 = primary.J2;
Re = primary.Re;
mu = primary.mu;
Tday = primary.Tday;
we = primary.we;
dOs = primary.dOs;

% Handle Input
rA = hA + primary.Re;

% N-R intialization
tol = 1e-12;
delX = inf(3,1);
fVal = ones(3,1);
iter = 0;
maxIter = 20;

initSma = (((k*Tday/j)/2/pi)^2*mu)^(1/3);
initEcc = rA/initSma - 1;
initInc = acos(-dOs/(3/2*J2*Re^2*sqrt(mu)*initSma^(-7/2)*(1-initEcc^2)^-2));
x = [initEcc;
    initSma;
    initInc];

while norm(delX,inf) > tol && iter < maxIter
    J2Effects = 3/4*J2*Re^2/(1 - x(1)^2)^2 * ...
        (2*j/k*cos(x(3)) - sqrt(1 - x(1)^2)*(3*cos(x(3))^2 - 1) - (5*cos(x(3))^2 - 1));
    % Apogee Altitude
    f1 = rA - x(2)*(1 + x(1)); 
    % RGT Condition
    f2 = x(2)^2 - j*we*(pi/180)/(k*sqrt(mu))*x(2)^(7/2) - J2Effects;
    % Sun-synch condition
    f3 = dOs + 3/2*J2*Re^2*sqrt(mu)*x(2)^(-7/2)*(1-x(1)^2)^-2*cos(x(3));
    
    % Partial derivatives
    df1de = -x(2);
    
    df1da = -(1 + x(1));
    
    df1di = 0;
    
    df2de = 3/4*J2*Re^2*...
        x(1)/(1 - x(1)^2)^3*...
        (-4*(5*cos(x(3))^2 - 1) + 8*j/k*cos(x(3)) - sqrt(1-x(1)^2)*...
        3*(3*cos(x(3))^2 - 1));
    
    df2da = 2*x(2) - 7*j*we*(pi/180)/(2*k*sqrt(mu))*x(2)^(5/2);
    
    df2di = -3/4*J2*Re^2/(1 - x(1)^2)^2 * ...
        (-2*j/k*sin(x(3)) + ...
        sqrt(1 - x(1)^2)*6*cos(x(3))*sin(x(3)) + ...
        10*cos(x(3))*sin(x(3)));
    
    df3de = 6*J2*Re^2*sqrt(mu)*x(2)^(-7/2)*x(1)*(1-x(1)^2)^-2*cos(x(3));
    
    df3da = -21/4*J2*Re^2*sqrt(mu)*x(2)^(-9/2)*(1-x(1)^2)^-2*cos(x(3));
    
    df3di = -3/2*J2*Re^2*sqrt(mu)*x(2)^(-7/2)*(1-x(1)^2)^-2*sin(x(3));
    
    % Jacobian matrix & iteration
    J = [df1de df1da df1di;
         df2de df2da df2di;
         df3de df3da df3di];
    fVal = [f1; f2; f3];
    delX = -J\fVal;
    
    x = x + delX;
    iter = iter + 1;
end
if x(1) > 0 && isreal(x)
    ecc = x(1);
    sma = x(2);
    inc = x(3)*180/pi;
else
    error('No viable solution found!')
end
end