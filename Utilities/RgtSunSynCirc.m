function [ sma, ecc, inc, fVal ] = RgtSunSynCirc( j, k )
%RgtSunSynElements calculates i, e & a for an RGT, sun-synchronous Orbit 
%with given apogee height
% Setting hA to 0 will result in a circular orbit
%   Uses multivariate Newton - Raphson method

%   ~~~~~~~~~~~~~  INPUTS  ~~~~~~~~~~~~~~
%   j   - number of repeats
%   k   - number of days to repeat


% Planet
primary = earth();
J2 = primary.J2;
Re = primary.Re;
mu = primary.mu;
Tday = primary.Tday;
we = primary.we;
dOs = primary.dOs;

% N-R intialization
tol = 1e-12;
delX = inf(2,1);
fVal = ones(2,1);
iter = 0;
maxIter = 20;

initSma = (((k*Tday/j)/2/pi)^2*mu)^(1/3);
initInc = acos(-dOs/(3/2*J2*Re^2*sqrt(mu)*initSma^(-7/2)));
x = [initSma;
    initInc];

while norm(delX,inf) > tol && iter < maxIter
    J2Effects = 3/4*J2*Re^2 * ...
        (2*j/k*cos(x(2)) - (3*cos(x(2))^2 - 1) - (5*cos(x(2))^2 - 1));

    % RGT Condition
    f1 = x(1)^2 - j*we*(pi/180)/(k*sqrt(mu))*x(1)^(7/2) - J2Effects;
    % Sun-synch condition
    f2 = dOs + 3/2*J2*Re^2*sqrt(mu)*x(1)^(-7/2)*cos(x(2));
    
    % Partial derivatives
    df1da = 2*x(1) - 7*j*we*(pi/180)/(2*k*sqrt(mu))*x(1)^(5/2);
    
    df1di = -3/4*J2*Re^2 * ...
        (-2*j/k*sin(x(2)) + ...
        6*cos(x(2))*sin(x(2)) + ...
        10*cos(x(2))*sin(x(2)));
    
    df2da = -21/4*J2*Re^2*sqrt(mu)*x(1)^(-9/2)*cos(x(2));
    
    df2di = -3/2*J2*Re^2*sqrt(mu)*x(1)^(-7/2)*sin(x(2));
    
    % Jacobian matrix & iteration
    J = [df1da df1di;
         df2da df2di];
    fVal = [f1; f2];
    delX = -J\fVal;
    
    x = x + delX;
    iter = iter + 1;
end
if x(1) > 0 && isreal(x)
    ecc = 0;
    sma = x(1);
    inc = x(2)*180/pi;
else
    error('No viable solution found!')
end
end