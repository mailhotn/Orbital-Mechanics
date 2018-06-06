function [ x ] = CalcRgtElement( a, e, inc, j, k, primary , tol)
%CalcRgtElement calculates the missing element for a RGT orbit
%   The missing element is left empty ([])  as input
%   ~~~~~~~~~~~~~  INPUTS  ~~~~~~~~~~~~~~
%   a     - semimajor axis [km]
%   e     - eccentricity
%   inc   - orbital inclination [deg]
%   j     - number of repeats
%   k     - number of days to repeat

if nargin < 6 % default is Earth
    primary = earth();
end
if nargin < 7 % default is Earth
    tol = 1e-12;
end
ratio = inf;
iter = 0;
maxIter = 20;
if isempty(a)
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
    
elseif isempty(e) % Numerical Issues
    e = 0.01;
    while abs(ratio) > tol && iter < maxIter % Newton - Raphson
        f  = primary.mu^(1/2)/a^(3/2) - (j*primary.we*(pi/180))/k + ...
            (3*primary.J2*primary.Re^2*(5*cosd(inc)^2 + (1 - e^2)^(1/2)*...
            (3*cosd(inc)^2 - 1) - (2*j*cosd(inc))/k - 1))/(4*(1 - e^2)^2);
        df = -(3*primary.J2*primary.Re^2*e*(3*cosd(inc)^2 - 1))/...
            (4*(1 - e^2)^(5/2)) - (3*primary.J2*primary.Re^2*e*(5*cosd(inc)^2 +...
            (1 - e^2)^(1/2)*(3*cosd(inc)^2 - 1) - (2*j*cosd(inc))/k - 1))/(1 - e^2)^3;
        ratio = f/df;
        e = e - ratio;
        iter = iter + 1;
    end
    if imag(e) == 0 && real(e) > 0
        x = e;
    else
        error(['No viable solution found!'])
    end
    
elseif isempty(inc) % Numerical Issues
    inc = pi/4;
    while abs(ratio) > tol && iter < maxIter % Newton - Raphson
        f  = primary.mu^(1/2)/a^(3/2) - (j*primary.we*(pi/180))/k + ...
            (3*primary.J2*primary.Re^2*(5*cos(inc)^2 + (1 - e^2)^(1/2)*...
            (3*cos(inc)^2 - 1) - (2*j*cos(inc))/k - 1))/(4*(e^2 - 1)^2);
        df = -(3*primary.J2*primary.Re^2*(10*cos(inc)*sin(inc) - ...
            (2*j*sin(inc))/k + 6*cos(inc)*sin(inc)*(1 - e^2)^(1/2)))/...
            (4*(e^2 - 1)^2);
        ratio = f/df;
        inc = inc - ratio;
        iter = iter + 1;
    end
    if imag(inc) == 0 && real(inc) > 0
        x = wrapTo180(inc*180/pi);
    else
        error(['No viable solution found!']) %#ok<*NBRAK>
    end
else
    error('Orbit is already defined')
end
end

