syms mu a tau wE inc Re J2 altp real

p = 2*(Re+altp)-(Re+altp)^2/a;
e = 1-(Re+altp)/a;
n = sqrt(mu/a^3);

xi = 3*Re^2*J2/4*p^-2;
chi = 1 + xi*(4 + 2*sqrt(1-e^2) - (5 + 3*sqrt(1-e^2))*sin(inc)^2);

f = mu^(-1/2)*a^(3/2) - tau/wE*chi/(1 + 2/wE*xi*n*cos(inc));

simplify(diff(f,a))

%%
Earth = earth();
tol = 1e-10;
iter = 0;
ratio = inf;

j = 51;
k = 10;
tau = k/j;
inc = 120*pi/180;
hp = 500;

mu = Earth.mu;
we = Earth.we*pi/180;
re = Earth.Re;
J2 = Earth.J2;

a = ((tau/we)^2*mu)^(1/3);

while abs(ratio) > tol
    p = 2*(re+hp)-(re+hp)^2/a;
    e = 1-(re+hp)/a;
    n = sqrt(mu/a^3);
    xi = 3*re^2*J2/4*p^-2;
    chi = 1 + xi*(4+2*sqrt(1-e^2)-(5+3*sqrt(1-e^2))*sin(inc)^2);
    
    f = mu^(-1/2)*a^(3/2) - tau/we*chi/(1 + 2/we*xi*n*cos(inc));
    df = (3*a^(1/2))/(2*mu^(1/2)) - (tau*((3*J2*re^2*(4*a + 2*sign(a)*(-(hp + re)*(hp - 2*a + re))^(1/2) - 5*a*sin(inc)^2 - 3*sign(a)*sin(inc)^2*(-(hp + re)*(hp - 2*a + re))^(1/2)))/(2*(hp + re)*(hp - 2*a + re)^3) - (3*J2*re^2*sign(a)*(3*sin(inc)^2 - 2)*(hp - a + re))/(4*(hp + re)*(-(hp + re)*(hp - 2*a + re))^(1/2)*(hp - 2*a + re)^2)))/(we*((3*J2*a^2*re^2*cos(inc)*(mu/a^3)^(1/2))/(2*we*(hp + re)^2*(hp - 2*a + re)^2) + 1)) + (3*J2*mu*re^2*tau*cos(inc)*((3*J2*a*re^2*(4*a + 2*sign(a)*(-(hp + re)*(hp - 2*a + re))^(1/2) - 5*a*sin(inc)^2 - 3*sign(a)*sin(inc)^2*(-(hp + re)*(hp - 2*a + re))^(1/2)))/(4*(hp + re)^2*(hp - 2*a + re)^2) + 1)*(6*a + hp + re))/(4*a^2*we^2*((3*J2*a^2*re^2*cos(inc)*(mu/a^3)^(1/2))/(2*we*(hp + re)^2*(hp - 2*a + re)^2) + 1)^2*(hp + re)^2*(mu/a^3)^(1/2)*(hp - 2*a + re)^3);
    ratio = f/df;
    a = a - ratio;
    iter = iter + 1;
end
e = 1 - (re+hp)/a;
a2 = CalcRgtElement([],e,inc*180/pi,j,k);