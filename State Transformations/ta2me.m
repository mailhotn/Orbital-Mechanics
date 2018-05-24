function [ Me ] = ta2me(th, e)
%true2mean_anomaly calculates the Mean anomaly of a spacecraft given it's
%true anomaly and eccentricity.  All angles are degrees.
th = wrapTo2Pi(th*pi/180);

E = 2*atan(sqrt((1-e)./(1+e)).*tan(th/2));
Me = E - e.*sin(E);

Me = wrapTo360(Me*180/pi);
end