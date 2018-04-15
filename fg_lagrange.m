function [ R, V ] = fg_lagrange( R0, V0, dt, tol, mu )
%fg_lagrange calculates change of position & velocity over time interval
%   calculates lagrange coefficients f, g, f_dot, g_dot and applies them to
%   the initial state to get the new state after dt
%   Solves for universal variable X using Newton-Raphson method.
%
%  ~~~~~~~~~~   Variables    ~~~~~~~~~~~
% mu - gravitational parameter (km^3/s^2)
% r0 - the initial position vector (km)
% v0 - the initial velocity vector (km/s)
% r  - the final position vector (km)
% v  - the final velocity vector (km/s)
% dt - elapsed time (s)
if nargin < 4
    tol = 1e-14;
end
if nargin < 5
    mu = 398600; % default is earth orbit of satellite with negligible mass
end
r0 = norm(R0);
v0 = norm(V0);
vr0 = dot(V0,R0)/r0;
a = 2/r0 - v0^2/mu; % inverse of semimajor axis
X = sqrt(mu)*abs(a)*dt;
ratio = inf;
while ratio > tol
    f = (r0*vr0/sqrt(mu))*X^2*StumC(a*X^2) +...
        (1-a*r0)*X^3*StumS(a*X^2) + r0*X - sqrt(mu)*dt;
    df = (r0*vr0/sqrt(mu))*X*(1 - a*X^2*StumS(a*X^2)) +...
        (1 - a*r0)*X^2*StumC(a*X^2) + r0;
    ratio = f/df;
    X = X - ratio;
end
f = 1 - X^2/r0*StumC(a*X^2);
g = dt - 1/sqrt(mu)*X^3*StumS(a*X^2);
R = f*R0 + g*V0;
rm = norm(R);
df = sqrt(mu)/(rm*r0)*(a*X^3*StumS(a*X^2)-X);
dg = 1 - X^2/rm*StumC(a*X^2);
V = df*R0 + dg*V0;
end