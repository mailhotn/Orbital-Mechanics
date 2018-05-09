function [ Me ] = ta2me(th, e)
%true2mean_anomaly calculates the Mean anomaly of a spacecraft given it's
%true anomaly and eccentricity.  All angles are degrees.
th = th*pi/180;

E = 2*atan(sqrt((1-e)./(1+e)).*tan(th/2));
Me = wrapTo2Pi(E - e.*sin(E));

Me = Me*180/pi;
end