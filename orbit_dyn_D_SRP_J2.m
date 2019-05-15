function [ dX ] = orbit_dyn_D_SRP_J2( t, X )
% orbit_dyn_D_SRP_J2 is used to simulate the dynamics of a spacecraft in
% earth orbit with Drag, SRP, and J2 perturbations.
% The function outputs the derivative of the state vector for use in a
% numerical integrator such as ode45
global K_d;
global K_s;
global day;
mu = 398600.440;
J2 = 0.0010826265;
Re = 6378;

R = X(1:3);
V = X(4:6);
r = norm(R);
v = norm(V);
f_d = -0.5*EarthAtmosphere(r)*v*K_d*V;
f_s = -4.55e-6*K_s*sun_dir(day);
OE = eci2oe(R,V);
inc = OE(3,:);
w   = OE(5,:);
th  = OE(6,:);
phi = asin(sind(th+w)*sind(inc));
f_J2 = -mu*J2*Re^2/r^4*(3*sin(phi)*[0 0 1].' + (-7.5*sin(phi)^2+1.5)*R./r);

dX = [V; -mu/r^3*R + f_d + f_s + f_J2];
end

