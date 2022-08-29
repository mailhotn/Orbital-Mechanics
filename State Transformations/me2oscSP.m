function [OE_osc] = me2oscSP(OE_m,primary)
%me2oscSP Calculates the osculating orbital elements of a satellite
%   Based on Kozai 1959
%   Accepts a 6xN Matrix of osculating elements in the following order:
%
% a  - semimajor axis (km)
% e  - eccentricity
% i  - inclination relative to equatorial plane (deg)
% O  - right ascension of ascending node (deg)
% w  - argument of periapsis (deg)
% Me - Mean anomaly (deg)
% Outputs a 6xN Matrix on mean elements in the same order

if nargin == 1
    primary = earth();
end
    

J2 = primary.J2;
Re = primary.Re;

sma0 = OE_m(1,:);
ecc0 = OE_m(2,:);
inc0 = pi/180*OE_m(3,:);
ran0 = wrapTo2Pi(pi/180*OE_m(4,:));
aop0 = wrapTo2Pi(pi/180*OE_m(5,:));
man0 = wrapTo360(OE_m(6,:));
f = me2ta(man0,ecc0);
man0 = pi/180*man0;
f = pi/180*f;

A2 = J2*Re^2*1.5; % positive A2: Me 2 Osc

dSmaS = 