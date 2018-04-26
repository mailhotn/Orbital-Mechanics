function [ Me ] = ta2me(th, e)
%true2mean_anomaly calculates the Mean anomaly of a spacecraft given it's
%true anomaly and eccentricity.  All angles are radians.
E = 2*atan(sqrt((1-e)./(1+e)).*tan(th/2));
Me = E - e.*sin(E);
if Me < 0
    Me = 2*pi + Me;
end
end