function [ r, v ] = fg_lagrange( r0, v0, dt, mu )
%fg_lagrange calculates change of position & velocity over time interval
%   calculates lagrange coefficients f, g, f_dot, g_dot and applies them to
%   the initial state to get the new state after dt
%
%  ~~~~~~~~~~   Variables    ~~~~~~~~~~~
% mu - gravitational parameter (km^3/s^2)
% r0 - the initial position vector (km)
% v0 - the initial velocity vector (km/s)
% r  - the final position vector (km)
% v  - the final velocity vector (km/s)
% dt - elapsed time (s)

if nargin < 4
    mu = 398600; % default is earth orbit of satellite with negligible mass
end
syms X;
r0m = norm(r0);
v0m = norm(v0);
vr0 = dot(v0,r0)/r0m;
a = 2/r0m - v0m^2/mu; % inverse of semimajor axis
eq = sqrt(mu)*dt == (r0m*vr0/sqrt(mu))*X^2*StumC(a*X^2,a) + (1-a*r0m)*X^3*StumS(a*X^2,a) + r0m*X; %#ok
X = double(vpasolve(eq,X,sqrt(mu)*abs(a)*dt));
f = 1 - X^2/r0m*StumC(a*X^2);
g = dt - 1/sqrt(mu)*X^3*StumS(a*X^2);
r = f*r0 + g*v0;
rm = norm(r);
df = sqrt(mu)/(rm*r0m)*(a*X^3*StumS(a*X^2)-X);
dg = 1 - X^2/rm*StumC(a*X^2);
v = df*r0 + dg*v0;
end